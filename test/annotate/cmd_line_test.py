'''
This script tests methods and classes in unifam.annotate.cmd_line module

Ideally we should do unit test, for now this is just some adhoc test.
'''
import os
import sys


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


from unifam.annotate.cmd_line import CmdLineGenerator
from unifam.base.util import SysUtil

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