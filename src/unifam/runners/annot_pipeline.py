'''
Annotatation pipeline.
This pipeline reads the output files from various programs, parse them,
generate input files for pathway tools, and run it.

Sample command
--------------
# suppose unifam repository is at ~/git/unifam

time PYTHONPATH=~/git/unifam/src python ~/git/unifam/src/unifam/runners/annot_pipeline.py \
        --config_file=~/git/unifam/examples/annotation/prok_isolate_genome.cfg

# modify the config for file locations before running
'''
import argparse
import sys
import logging
import os
from configparser import ConfigParser
from unifam.annotate.cmd_line import CmdLineGenerator
from unifam.annotate.parse_util import ProteinAnnotator
from unifam.annotate.parse_util import ProdigalOutputParser
from unifam.annotate.dir_spec import UniFamDirSpec
from unifam.base.util import SysUtil
from unifam.db.swiss_prot_parser import SwissProtParser
from unifam.annotate.parse_util import RnaOutputReader
from unifam.annotate.pathologic_input import PathoLogicInput


def create_parser():
    parser = argparse.ArgumentParser(description="Annotate usearch protein clusters",
                                     prog='annot_cluster',
                                     prefix_chars='-',
                                     fromfile_prefix_chars='@',
                                     conflict_handler='resolve',
                                     add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--config_file", required=True,
                        help="config file including all the input and paramter specifications")
    return parser


def set_up_logger():
    # TODO: set up file logger with specified log file name, and blackhole logger

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    #file_handler = logging.FileHandler('unifam.log')
    #file_handler.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(filename)s - %(levelname)s - %(lineno)d'
                                  '%(process)d: %(message)s')
    #file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    #logger.addHandler(file_handler)
    logger.addHandler(console_handler)


def main():

    set_up_logger()
    parser = create_parser()
    args = parser.parse_args()
    config_file = args.config_file
    assert os.path.isfile(config_file), f'{config_file} does not exist'

    config = ConfigParser()
    config.read(config_file)
    assert isinstance(config, ConfigParser), type(config)
    unifam_dir_spec = UniFamDirSpec.from_config(config)

    cmd_line_gen = CmdLineGenerator(config_file)

    prog_to_cmd = {'prodigal': cmd_line_gen.get_prodigal_cmd(),
                   'rnammer': cmd_line_gen.get_RNAmmer_cmd(),
                   'tRNAscan': cmd_line_gen.get_tRNAscan_cmd()
                   }

    for prog, prog_cmd in prog_to_cmd.items():
        print(f'>>> Run {prog}:')
        logging.info(f'>>> Run {prog}:')
        SysUtil.run_cmd(prog_cmd, timeout=None, dry_run=False)

    # hmmsearch uses prodigal output faa file, so we run it separately
    # after protein, tRNA, and rRNA are found
    print('>>> Run hmmsearch')
    logging.info('>>> Run hmmsearch')
    hmmsearch_cmd = cmd_line_gen.get_hmmsearch_cmd()
    SysUtil.run_cmd(hmmsearch_cmd, timeout=None, dry_run=False)

    # load cluster annotation
    cluster_annot_path = os.path.expandvars(config.get('UniFam', 'cluster_annot_file'))
    logging.info(f'Load cluster annotation file to dict from {cluster_annot_path}')
    annot_dict = SwissProtParser.load_annot(cluster_annot_path)
    logging.info(f'Cluster annot dict size: {len(annot_dict)}')

    # hmmsearch result
    domtab_file = unifam_dir_spec.get_hmmsearch_domtbl_file()
    seq_coverage = config.getfloat('UniFam', 'seq_coverage', fallback=0.5)
    hmm_coverage = config.getfloat('UniFam', 'hmm_coverage', fallback=0.5)
    domtab_df = ProteinAnnotator.parse_domtbl_file(domtab_file, seq_coverage, hmm_coverage)
    logging.info(f'Parsed domtab output file {domtab_file} to data frame with length {len(domtab_df)}')

    # contig information from prodigal .gbk output
    gbk_file = unifam_dir_spec.get_prodigal_out_file(out_type='gbk')
    contig_to_info = ProdigalOutputParser.get_contig_info(gbk_file)
    logging.info(f'Parsed gbk file {gbk_file} to {len(contig_to_info)} contigs')

    # tRNA and rRNA results
    rrna_gff_file = unifam_dir_spec.get_RNAmmer_gff_file()
    rrna_df = RnaOutputReader.read_RNAmmer_gff(rrna_gff_file)
    logging.info(f'Parsed rRNA gff file {rrna_gff_file} to data frame with length {len(rrna_df)}')
    trna_o_file = unifam_dir_spec.get_tRNAscan_output_file()
    trna_df = RnaOutputReader.read_tRNAscan_output(trna_o_file)
    logging.info(f'Parsed tRNAscan .o file {trna_o_file} to data frame with length {len(trna_df)}')

    # write all input files for pathologic
    patho_input_dir = unifam_dir_spec.get_patho_input_dir()
    SysUtil.mkdir_p(patho_input_dir)
    genetic_element_file = os.path.join(patho_input_dir, 'genetic-elements.dat')
    PathoLogicInput.write_genetic_element(genetic_element_file, contig_to_info, domtab_df, annot_dict, rrna_df, trna_df)

    # We use placeholder for the dbID and name in the example below
    organism = config.get('PathoLogic', 'organism')
    domain = config.get('PathoLogic', 'domain')
    db_id = config.get('PathoLogic', 'db_id')
    codon_table = config.get('PathoLogic', 'codon_table', fallback='11')
    PathoLogicInput.write_organism_params(os.path.join(patho_input_dir, 'organism-params.dat'),
                                          PathoLogicInput.get_organism_params_dict(organism, db_id, domain=domain,
                                                                                   codon_table=codon_table))
    logging.info(f'Written pathway-tools input files to {patho_input_dir}')

    pathologic_cmd = cmd_line_gen.get_pathologic_cmd()
    print(f'>>> Run pathway-tools:')
    SysUtil.run_cmd(pathologic_cmd, timeout=None, dry_run=False)


if __name__ == '__main__':
    sys.exit(main())
