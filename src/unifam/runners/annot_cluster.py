'''
Annotate protein clusters found by usearch,
and save the annotation to disk.

Sample command
--------------
# suppose unifam repository is at ~/git/unifam

time PYTHONPATH=~/git/unifam/src python ~/git/unifam/src/unifam/runners/annot_cluster.py \
        --sp_json_file=uniprot_sprot.json --uc_file=uniprot_nomsa.uc --out_file=cluster_annot.json \
        --share_thresh=0.5

# time taken on a desktop machine
...
...
annotated 270000 clusters
Save cluster annotation to json format at cluster_annot.json

real    289m36.492s
user    289m27.548s
sys     0m8.416s

'''
import argparse
import sys
import json
from unifam.db.usearch_helper import UsearchCluster
from unifam.db.swiss_prot_parser import SwissProtParser
from unifam.db.swiss_prot_parser import ClusterAnnot


def create_parser():
    parser = argparse.ArgumentParser(description="Annotate usearch protein clusters",
                                     prog='annot_cluster',
                                     prefix_chars='-',
                                     fromfile_prefix_chars='@',
                                     conflict_handler='resolve',
                                     add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sp_json_file", required=True, help="swissprot protein annotation json file")
    parser.add_argument("--uc_file", required=True, help="Usearch cluster file")
    parser.add_argument("--out_file", required=True, help="path to save the cluster annot dict")
    parser.add_argument("--share_thresh", type=float,
                        help="real number in [0,1], threshold used to aggregate protein annotation")
    return parser


def main():

    parser = create_parser()
    args = parser.parse_args()
    print(f'Load swissprot annotation from json file: {args.sp_json_file}')
    id_to_annot_dict = SwissProtParser.load_annot(args.sp_json_file)
    print(f'Load usearch cluster file: {args.uc_file}')
    uc_parser = UsearchCluster(args.uc_file)
    cluster_to_annot = dict()
    print(f'Annotate protein clusters')
    for idx, cluster in enumerate(sorted(list(uc_parser.get_sp_cluster_set()))):
        if idx > 0 and idx % 1000 == 0:
            print(f'annotated {idx} clusters')
        seq_id_list = uc_parser.get_seqs_in_cluster(cluster)
        cluster_annot = ClusterAnnot([id_to_annot_dict[seq_id] for seq_id in seq_id_list])
        cluster_to_annot[cluster] = cluster_annot.get_cluster_annot_dict(args.share_thresh)
    print(f'Save cluster annotation to json format at {args.out_file}')
    with open(args.out_file, 'w') as save_f:
        json.dump(cluster_to_annot, save_f, indent=2)


if __name__ == '__main__':
    sys.exit(main())
