'''
Parse swissprot .dat file, extract part of the annotation and
save the resulted dict to a json formatted file

Sample command
--------------
# suppose unifam repository is at ~/git/unifam

time PYTHONPATH=~/git/unifam/src python ~/git/unifam/src/unifam/runners/parse_swissprot_annot.py \
        --dat_file=uniprot_sprot.dat --out_file=uniprot_sprot.json
'''
import argparse
import sys
import os
from unifam.db.swiss_prot_parser import SwissProtParser


def create_parser():
    parser = argparse.ArgumentParser(description="Parse swiss prot .dat file and save to json",
                                     prog='parse_swissprot_annot',
                                     prefix_chars='-',
                                     fromfile_prefix_chars='@',
                                     conflict_handler='resolve',
                                     add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--dat_file", required=True, help="swissprot .dat file path")
    parser.add_argument("--out_file", required=True, help="path to save the parsed annot dict")
    return parser


def main():

    parser = create_parser()
    args = parser.parse_args()
    assert not os.path.isdir(args.out_file), f'{args.out_file} is an existing dir'
    id_to_annot = SwissProtParser.read_annot_file(args.dat_file)
    SwissProtParser.save_annot(id_to_annot, args.out_file)


if __name__ == '__main__':
    sys.exit(main())
