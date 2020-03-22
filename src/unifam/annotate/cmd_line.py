'''
command line related classes and functions
'''
import ast
from configparser import ConfigParser
import os

from unifam.annotate.dir_spec import UniFamDirSpec
from unifam.base.util import SysUtil


class CmdLineGenerator(object):
    """
    Class to generate command lines from config file, using ConfigParser

    Parameters
    ----------
    config_file
        config file that can be read by ConfigParser
    """

    def __init__(self, config_file):
        assert os.path.isfile(config_file), config_file
        self._config = ConfigParser()
        self._config.read(config_file)
        self._unifam_dir_spec = UniFamDirSpec.from_config(self._config)
        self._input_file = self._config.get('UniFam', 'input_file')

    def get_prodigal_cmd(self):
        assert self._config.has_section('prodigal')
        prodigal_dict = dict(self._config['prodigal'])
        prodigal_dict.update({'input_file': self._input_file,
                              'translation_file': self._unifam_dir_spec.get_prodigal_out_file('faa'),
                              'output_file': self._unifam_dir_spec.get_prodigal_out_file('gbk')})
        boolean_key_list = ['allow_run_off_edge', 'mask_n_runs', 'quiet']
        for bkey in boolean_key_list:
            prodigal_dict[bkey] = ast.literal_eval(prodigal_dict[bkey])
        return CmdLineHelper.get_prodigal_cmd(**prodigal_dict)

    def get_RNAmmer_cmd(self):
        assert self._config.has_section('RNAmmer')
        RNAmmer_dict = dict(self._config['RNAmmer'])
        RNAmmer_dict.update({'input_file': self._input_file, 'domain': self._config.get('UniFam', 'domain'),
                             'output_gff_file': self._unifam_dir_spec.get_RNAmmer_gff_file()})
        return CmdLineHelper.get_RNAmmer_cmd(**RNAmmer_dict)

    def get_tRNAscan_cmd(self):
        assert self._config.has_section('tRNAscan')
        tRNAscan_dict = dict(self._config['tRNAscan'])
        if 'quiet' in tRNAscan_dict:
            tRNAscan_dict['quiet'] = ast.literal_eval(tRNAscan_dict['quiet'])
        tRNAscan_dict.update({'input_file': self._input_file, 'domain': self._config.get('UniFam', 'domain'),
                              'out_prefix': self._unifam_dir_spec.get_prefix()})
        return CmdLineHelper.get_tRNAscan_cmd(**tRNAscan_dict)

    def get_hmmsearch_cmd(self):
        assert self._config.has_section('hmmsearch')
        hmmsearch_dict = dict(self._config['hmmsearch'])
        hmmsearch_dict.update({'input_file': self._input_file,
                               'domtbl_out_file': self._unifam_dir_spec.get_hmmsearch_domtbl_file(),
                               'output_file': self._unifam_dir_spec.get_hmmsearch_out_file()})
        return CmdLineHelper.get_hmmsearch_cmd(**hmmsearch_dict)


class CmdLineHelper(object):
    """
    Class used to help assemble command lines to run in shell, and execute them.
    """
    @classmethod
    def get_prodigal_cmd(cls, prodigal_path, input_file, translation_file, output_file,
                         allow_run_off_edge=True, nucleotide_seq_file=None, translation_table='11',
                         mask_n_runs=False, procedure='single', quiet=False, protein_n_score_file=None):
        """
        Returns the command line to run prodigal from given specification, and make parent directoy for the output files.

        prodigal -h

        Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
                         [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
                         [-p mode] [-q] [-s start_file] [-t training_file] [-v]

         -a:  Write protein translations to the selected file.
         -c:  Closed ends.  Do not allow genes to run off edges.
         -d:  Write nucleotide sequences of genes to the selected file.
         -f:  Select output format (gbk, gff, or sco).  Default is gbk.
         -g:  Specify a translation table to use (default 11).
         -h:  Print help menu and exit.
         -i:  Specify FASTA/Genbank input file (default reads from stdin).
         -m:  Treat runs of N as masked sequence; don't build genes across them.
         -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
         -o:  Specify output file (default writes to stdout).
         -p:  Select procedure (single or meta).  Default is single.
         -q:  Run quietly (suppress normal stderr output).
         -s:  Write all potential genes (with scores) to the selected file.
         -t:  Write a training file (if none exists); otherwise, read and use
              the specified training file.
         -v:  Print version number and exit.


        """
        assert SysUtil.is_executable_file(prodigal_path), f'{prodigal_path} is not an executable file'
        assert os.path.isfile(input_file), f'{input_file} does not exist'
        SysUtil.mkdir_p(os.path.dirname(translation_file))
        SysUtil.mkdir_p(os.path.dirname(output_file))
        assert isinstance(allow_run_off_edge, bool)
        close_ends_tag = '' if allow_run_off_edge else '-c'
        # -d
        nucleotide_file_tag = ''
        if nucleotide_seq_file is not None:
            SysUtil.mkdir_p(nucleotide_seq_file)
            nucleotide_file_tag = '-d {nucleotide_seq_file}'

        mask_N_runs_tag = '-m' if mask_n_runs else ''
        assert procedure in {'single', 'meta'}, procedure
        quiet_tag = '' if not quiet else '-q'
        protein_and_score_file_tag = ''
        if protein_n_score_file is not None:
            SysUtil.mkdir_p(os.path.dirname(protein_n_score_file))
            protein_and_score_file_tag = '-s {protein_n_score_file}'

        return (f'{prodigal_path} -i {input_file} -a {translation_file} -o {output_file} '
                f'{nucleotide_file_tag} -f gbk -g {translation_table} {mask_N_runs_tag} {close_ends_tag} '
                f'-p {procedure} {quiet_tag} {protein_and_score_file_tag}')

    @classmethod
    def get_RNAmmer_cmd(cls, rnammer_path, input_file, domain, output_gff_file):
        """
        Returns the command line to run rnammer from given specifications, and make parent directoy for the output files.

        [rnammer man page](https://www.cbs.dtu.dk/cgi-bin/nph-runsafe?man=rnammer)

        RNAmmer - predicts ribosomal RNA genes in prokaryotic genome sequences

        SYNOPSIS
               rnammer [-S kingdom] [-m molecules]  [-xml xml-file] [-gff gff-file] [-d] [-p]
               [-h hmmreport] [-f fasta-file] [sequence]

               or

               rnammer [-S kingdom] [-m molecules]  [-xml xml-file] [-gff gff-file] [-d] [-p]
               [-h hmmreport] [-f fasta-file] < [sequence]

        CORE
               core-rnammer [configuration]


        REQUIREMENTS
               The main executable 'rnammer' requires the core RNAmmer program 'core-rnammer'. The core program requires the
               binary 'hmmsearch'  (http://hmmer.org/) version 2.3.2

        OPTIONS
               -S     Specifies the super kingdom of the input sequence. Can be either 'arc', 'bac', or 'euk'.

               -gff output gff file
                      Specifies filename for output in GFF version 2 output

               -multi Runs all molecules and both strands in parallel

               -f fasta
                      Specifies filename for output fasta file of predicted rRNA genes

               -h hmmreport
                      Specifies filename for output HMM report.

               -m     Molecule type can be 'tsu' for 5/8s rRNA, 'ssu' for 16/18s rRNA, 'lsu' for 23/28s rRNA or
                      any combination separated by comma.

               [sequence]
                      The input file to process.

        """
        assert SysUtil.is_executable_file(rnammer_path), f'{rnammer_path} is not an executable file'
        assert os.path.isfile(input_file), input_file
        assert domain in {'arc', 'bac', 'euk'}, domain
        SysUtil.mkdir_p(os.path.dirname(output_gff_file))
        # Here we only need the gff file, so other options are not coded.
        return f'{rnammer_path} -S {domain} -m ssu,lsu,tsu -gff {output_gff_file} {input_file}'

    @classmethod
    def get_tRNAscan_cmd(cls, trnascan_path, input_file, domain, out_prefix, quiet=False):
        """
        Returns the command line to run tRNAscan from given specification, and make parent directoy for the output files.

        tRNAscan-SE-2.0 $ ./tRNAscan-SE -h
        Usage: tRNAscan-SE [-options] <FASTA file(s)>

          Scan a sequence file for tRNAs
           -- default: use Infernal & tRNA covariance models
              with eukaryotic sequences
              (use 'Search Mode Options' below to scan other types of sequences)

        Search Mode Options:

          -E                          : search for eukaryotic tRNAs (default)
          -B                          : search for bacterial tRNAs
          -A                          : search for archaeal tRNAs
          -M <model>                  : search for mitochondrial tRNAs
                                          options: mammal, vert
          -O                          : search for other organellar tRNAs
          -G                          : use general tRNA model (cytoslic tRNAs from all 3 domains included)
          --mt <model>                : use mito tRNA models for cytosolic/mito detemination
                                          (if not specified, only cytosolic isotype-specific model scan will be performed)
          -I                          : search using Infernal
                                          default use with -E, -B, -A, or -G; optional for -O
              --max                   : maximum sensitivity mode - search using Infernal without hmm filter (very slow)
          -L                          : search using the legacy method (tRNAscan, EufindtRNA, and COVE)
                                          use with -E, -B, -A or -G
          -C  --cove                  : search using COVE analysis only (legacy, extremely slow)
                                          default use with -O
          -H  --breakdown             : show breakdown of primary and secondary structure components to
                                          covariance model bit scores
          -D  --nopseudo              : disable pseudogene checking

        Output options:

          -o  --output <file>         : save final results in <file>
          -f  --struct <file>         : save tRNA secondary structures to <file>
          -s  --isospecific <file>    : save results using isotype-specific models in <file>
          -m  --stats <file>          : save statistics summary for run in <file>
                                          (speed, # tRNAs found in each part of search, etc)
          -b  --bed <file>            : save results in BED file format of <file>
          -a  --fasta <file>          : save predicted tRNA sequences in FASTA file format of <file>
          -l  --log <file>            : save log of program progress in <file>
          --detail                    : display prediction outputs in detailed view
          --brief                     : brief output format (no column headers)

          -? #                       : '#' in place of <file> chooses default name for output files
          -p  --prefix <label>        : use <label> prefix for all default output file names

          -d  --progress              : display program progress messages
          -q  --quiet                 : quiet mode (credits & run option selections suppressed)
          -y  --hitsrc                : show origin of hits (Ts=tRNAscan 1.4, Eu=EufindtRNA,
                                          Bo=Both Ts and Eu, Inf=Infernal)
        """
        assert SysUtil.is_executable_file(trnascan_path), f'{trnascan_path} is not an executable file'
        assert os.path.isfile(input_file), f'{input_file} is not an existing file'
        domain_to_option = {'bac': '-B', 'arc': '-A', 'euk': '-E'}
        domain_option = domain_to_option[domain]
        quiet_tag = '' if not quiet else '-q'
        assert isinstance(out_prefix, str), out_prefix
        SysUtil.mkdir_p(os.path.dirname(out_prefix))
        output_file = UniFamDirSpec.get_tRNAscan_out_file(out_prefix, 'o')
        struct_file = UniFamDirSpec.get_tRNAscan_out_file(out_prefix, 'structure')
        stat_file = UniFamDirSpec.get_tRNAscan_out_file(out_prefix, 'stat')

        return (f'{trnascan_path} {domain_option} '
                f'-o {output_file} -f {struct_file} -m {stat_file} {quiet_tag} '
                f'{input_file}')

    @classmethod
    def get_hmmsearch_cmd(cls, hmmsearch_path, input_file, database_dir, database_name, e_val_ub,
                          domtbl_out_file, output_file, n_cpu=1):
        """
        Returns the command line to run hmmsearch from given specification, and make parent directoy for the output files.

        Parameters
        ----------
        n_cpu : int
            number of cpu/threads to use, for option --cpu. Set to -1 to use default, i.e. not using this option

        $ ./hmmsearch -h
        # hmmsearch :: search profile(s) against a sequence database
        # HMMER 3.3 (Nov 2019); http://hmmer.org/
        # Copyright (C) 2019 Howard Hughes Medical Institute.
        # Freely distributed under the BSD open source license.
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Usage: hmmsearch [options] <hmmfile> <seqdb>

        Basic options:
          -h : show brief help on version and usage

        Options directing output:
          -o <f>           : direct output to file <f>, not stdout
          -A <f>           : save multiple alignment of all hits to file <f>
          --tblout <f>     : save parseable table of per-sequence hits to file <f>
          --domtblout <f>  : save parseable table of per-domain hits to file <f>
          --pfamtblout <f> : save table of hits and domains to file, in Pfam format <f>
          --acc            : prefer accessions over names in output
          --noali          : don't output alignments, so output is smaller
          --notextw        : unlimit ASCII text output line width
          --textw <n>      : set max width of ASCII text output lines  [120]  (n>=120)

        Options controlling reporting thresholds:
          -E <x>     : report sequences <= this E-value threshold in output  [10.0]  (x>0)
          -T <x>     : report sequences >= this score threshold in output
          --domE <x> : report domains <= this E-value threshold in output  [10.0]  (x>0)
          --domT <x> : report domains >= this score cutoff in output

        Options controlling inclusion (significance) thresholds:
          --incE <x>    : consider sequences <= this E-value threshold as significant
          --incT <x>    : consider sequences >= this score threshold as significant
          --incdomE <x> : consider domains <= this E-value threshold as significant
          --incdomT <x> : consider domains >= this score threshold as significant

        Options controlling model-specific thresholding:
          --cut_ga : use profile's GA gathering cutoffs to set all thresholding
          --cut_nc : use profile's NC noise cutoffs to set all thresholding
          --cut_tc : use profile's TC trusted cutoffs to set all thresholding

        Options controlling acceleration heuristics:
          --max    : Turn all heuristic filters off (less speed, more power)
          --F1 <x> : Stage 1 (MSV) threshold: promote hits w/ P <= F1  [0.02]
          --F2 <x> : Stage 2 (Vit) threshold: promote hits w/ P <= F2  [1e-3]
          --F3 <x> : Stage 3 (Fwd) threshold: promote hits w/ P <= F3  [1e-5]
          --nobias : turn off composition bias filter

        Other expert options:
          --nonull2     : turn off biased composition score corrections
          -Z <x>        : set # of comparisons done, for E-value calculation
          --domZ <x>    : set # of significant seqs, for domain E-value calculation
          --seed <n>    : set RNG seed to <n> (if 0: one-time arbitrary seed)  [42]
          --tformat <s> : assert target <seqfile> is in format <s>: no autodetection
          --cpu <n>     : number of parallel CPU workers to use for multithreads  [2]
        """
        assert SysUtil.is_executable_file(hmmsearch_path), f'{hmmsearch_path} is not an executable file'
        assert os.path.isfile(input_file), f'{input_file} is not an existing file'
        hmm_db_file = os.path.join(database_dir, f'{database_name}.hmm')
        assert os.path.isfile(hmm_db_file), hmm_db_file
        assert isinstance(n_cpu, int)
        SysUtil.mkdir_p(os.path.dirname(domtbl_out_file))
        SysUtil.mkdir_p(os.path.dirname(output_file))
        cpu_option = '' if n_cpu == -1 else f'--cpu {n_cpu}'

        return (f'{hmmsearch_path} -E {e_val_ub} --noali {cpu_option} --domtblout {domtbl_out_file} '
                f'-o {output_file} {input_file} {hmm_db_file}')
