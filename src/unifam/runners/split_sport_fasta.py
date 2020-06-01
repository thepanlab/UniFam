'''
Split uniprot SwissProt fasta file into 2 fasta files:
    1 with Eukaryote sequences
    1 with non-Eukaryote sequences

The domain is determined from json file including protein annotations parsed from uniprot_sprot.dat file.

Sample command
--------------
# suppose unifam repository is at ~/git/unifam

time PYTHONPATH=~/git/unifam/src python ~/git/unifam/src/unifam/runners/split_sport_fasta.py \
        --sp_json_file=uniprot_sprot.json --sp_fasta_file=uniprot_sprot.fasta --out_dir='.'

real    0m18.994s
user    0m13.064s
sys     0m2.968s
'''
import argparse
import sys
from unifam.db.swiss_prot_parser import SpSplitter


def create_parser():
    parser = argparse.ArgumentParser(description="Split SwissProt sequence file into 2 based on domain",
                                     prog='split_sprot_fasta',
                                     prefix_chars='-',
                                     fromfile_prefix_chars='@',
                                     conflict_handler='resolve',
                                     add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sp_json_file", required=True, help="Swissprot protein annotation json file")
    parser.add_argument("--sp_fasta_file", required=True, help="Swissprot protein fasta file")
    parser.add_argument("--out_dir", required=True, help="directory to save the 2 sub fasta files")
    return parser


def main():

    parser = create_parser()
    args = parser.parse_args()
    splitter = SpSplitter(args.sp_fasta_file, args.sp_json_file)
    splitter.split(args.out_dir)


if __name__ == '__main__':
    sys.exit(main())
