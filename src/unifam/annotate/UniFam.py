'''
Entrance point for running unifam pipeline.
'''
# Import Python modules
import argparse
import configparser
from datetime import datetime
import sys

from unifam.annotate.UniFam_lib import UniFam as unifam
from unifam.annotate.UniFam_lib_batch import UniFam as batch_unifam


# Import local modules
# Version
version_str = "1.1.0"
''' 0.0.1.  first version of pipeline, including prodigal, hmmsearch, and annotation
            all configuration options are in a file. User can(and must) change the options
            in the config file to customize.
    0.0.2.  Added more information for pathologic module of the analysis, and rRNA, tRNA
            analysis.
    1.0.0   First Stable version for release, UniFam 1.0.0
    1.1.0   Added README file to describe output files; zip pathway inference output results
    1.1.1   Working on migrating to python 3.6, and clean up code.
'''

parser = argparse.ArgumentParser(description="Annotation of contigs/proteins using UniFam",
                                 prog='UniFam',
                                 prefix_chars='-',
                                 fromfile_prefix_chars='@',
                                 conflict_handler='resolve',
                                 add_help=True,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                 )
# version control
parser.add_argument("--version", action="version", version='%(prog)s {}'.format(version_str))

# --verbose mode, default: false (quiet mode)
parser.add_argument("-v", "--verbose", action="store_true", help="verbose mode, more output")
parser.add_argument("-n", "--dryrun", action="store_true", help="dryrun, only print commands, do not execute")
parser.add_argument("-b", "--batch", action="store_true", help="batch mode to construct pathway", dest="batch")

# configuration file, required
parser.add_argument("-c", help="configuration file", dest='configFile', required=True)
# input file, required
parser.add_argument("-i", help="input fasta file (contig or protein fasta/faa file)", dest='inputFile', required=True)


def main(argv=None):

    if argv is None:
        args = parser.parse_args()

    # print some information
    if args.verbose:
        sys.stdout.write('Running verbosely\n'
                         f'configuration file = {args.configFile}\n'
                         f'input fasta file = {args.inputFile}\n'
                         f'batch mode = {args.batch}')

    # display work start, and time record
    start_time = datetime.now()
    sys.stderr.write("\n===============================================================================\n")
    sys.stderr.write(f"Welcome to UniFam v{version_str}: \n")

    # read configuration file
    config = configparser.ConfigParser()
    config.read(args.configFile)

    # Annotating with UniFam
    if args.batch:
        batch_unifam(args.inputFile, config, args.verbose)
    else:
        unifam(args.inputFile, config, args.verbose)

    # write the configuration file to standard output for checking
    # config.write(sys.stdout)

    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write("\nTotal Elapsed Time = [%s] \n" % duration)
    sys.stderr.write("===============================================================================\n")


if __name__ == '__main__':
    sys.exit(main())
