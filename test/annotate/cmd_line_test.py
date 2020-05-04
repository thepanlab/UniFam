'''
This script tests methods and classes in unifam.annotate.cmd_line module

Ideally we should do unit test, for now this is just some adhoc test.
'''
import logging
import os
import sys
import pandas as pd

from unifam.annotate.cmd_line import CmdLineGenerator
from unifam.annotate.parse_util import ProteinAnnotator
from unifam.base.util import SysUtil
from unifam.db.usearch_helper import UsearchCluster
from unifam.db.swiss_prot_parser import SwissProtParser


# First load python 3.6.6 on oscer
# module load Python/3.6.6-foss-2018b
# jcross@schooner1:unifam $ python --version
# Python 3.6.6
# Suppose you are at unifam root dir, Run this script with
# python test/annotate/cmd_line_test.py
# add {unifam}/src to sys.path
# this is not the best way to handle python path, we'll deal with this later
curr_file_path = __file__
unifam_root_dir = os.path.abspath(
    f'{os.path.dirname(os.path.dirname(os.path.dirname(curr_file_path)))}')
print(f'unifam repository root dir is {unifam_root_dir}')
sys.path.append(os.path.join(unifam_root_dir, 'src'))
print(f'Added  {sys.path[-1]} to sys.path')


# set up logger for logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('unifam.log')
file_handler.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(levelname)s - %(lineno)d'
                              '%(process)d: %(message)s')
file_handler.setFormatter(formatter)
console_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.addHandler(console_handler)

# find sample config file, and change working directory to its containing dir
prok_iso_genome_config_file = os.path.join(
    unifam_root_dir, 'examples', 'annotation', 'prok_isolate_genome.cfg')
os.chdir(os.path.dirname(prok_iso_genome_config_file))
print(f'current working directory is {os.getcwd()}')


cmd_line_gen = CmdLineGenerator(prok_iso_genome_config_file)
prodigal_cmd = cmd_line_gen.get_prodigal_cmd()
rnammer_cmd = cmd_line_gen.get_RNAmmer_cmd()
trnascan_cmd = cmd_line_gen.get_tRNAscan_cmd()
hmmsearch_cmd = cmd_line_gen.get_hmmsearch_cmd()

prog_to_cmd = {'prodigal': cmd_line_gen.get_prodigal_cmd(),
               'rnammer': cmd_line_gen.get_RNAmmer_cmd(),
               'tRNAscan': cmd_line_gen.get_tRNAscan_cmd(),
               'hmmsearch': cmd_line_gen.get_hmmsearch_cmd()
               }

# After inspecting all the commands and satisfied,
# change dry_run=True to False to actually run the commands
for prog, prog_cmd in prog_to_cmd.items():
    print(f'>>> Run {prog}:')
    SysUtil.run_cmd(prog_cmd, timeout=None, dry_run=True)

# ======================================================================

# Parsing Test
domtab_file = f'{unifam_root_dir}/example/prok_isolate_genome/output/GuestGenome.domtab'
full_df = ProteinAnnotator.read_domtab_file_to_df(domtab_file)  # all the records
df = ProteinAnnotator.parse_domtbl_file(domtab_file)  # Only the ones that pass the threshold, and choose the best hit
# parse result from old code
parse_out_file = f'{unifam_root_dir}/example/prok_isolate_genome/output/GuestGenome.domtab_parse.out'
out_df = pd.read_csv(parse_out_file, sep='\t', header=None, names=['seq_name', 'hmm_name', 'eval'])

# ======================================================================
# usearch cluster result parsing
uc_file = '/work/omicsbio/lizhang12/database_build/uniprot_nomsa.uc'
uc_parser = UsearchCluster(uc_file)
swiss_prot_cluster_list = sorted(list(uc_parser.get_sp_cluster_set()))
print(f'number of clusters containing at least 1 swiss prot sequence is {len(swiss_prot_cluster_list)}')
seq_list = uc_parser.get_seqs_in_cluster(swiss_prot_cluster_list[0])

# ======================================================================
# SwissProt .dat file parsing
dat_file = '/work/omicsbio/lizhang12/database_build/download/uniprot_sprot.dat'
id_to_annot = SwissProtParser.read_annot_file(dat_file, max_records=5000)
# Read 1000 records...
# Read 2000 records...
# Read 3000 records...
# Read 4000 records...
# Reached max: 5000 records, last read position: 33543260
# Read 5000 records, last read position: 33543260
