'''
Utility classes for parsing output
'''
import logging
import os

import numpy as np
import pandas as pd
from unifam.base.util import SysUtil


class ProteinAnnotator(object):
    """
    Class to annotate proteins found by hmmsearch with UniFam data base.
    """

    @classmethod
    def parse_domtbl_file(cls, domtbl_file, seq_coverage=0.5, hmm_coverage=0.5, output_file=None):
        """
        Parse the domain-wise table output from hmmsearch,
        find the best-matching models for all proteins in a genome.

        Things to consider:
        1) E-value
        2) Length of aligned region / total length of sequence and model > threshold (0.5).
        3) For each query sequence, return the record with the smallest E-value among all
           that pass the coverage threshold

        One can consider tighten the E-value threshold when running hmmsearch,
        to potentially speed up the search.

        Parameters
        ----------
        domtbl_file : str
            path to the domtbl output file from hmmsearch
            This file is fixed width format, with comment line starting with #
        seq_coverage : real
            lower limit of coverage of the query sequence
        hmm_coverage : real
            lower limit of coverage of the hmm
        output_file : str, optional
            path to the file where the domtab_df will be saved in csv format

        Returns
        -------
        pd.DataFrame
            data frame including the information on the best match of hmm from the input `domtbl_file`
            satisfying the specified conditions.
            Some important columns:
            seq_name, hmm_name, eval
        """
        assert seq_coverage > 0. and seq_coverage <= 1.0, seq_coverage
        assert hmm_coverage > 0. and hmm_coverage <= 1.0, hmm_coverage
        domtab_df = cls.read_domtab_file_to_df(domtbl_file).copy()
        min_coverage = min(seq_coverage, hmm_coverage)
        logging.info(f'Using seq_coverage >= {seq_coverage:.2f}, '
                     f'hmm_coverage >= {hmm_coverage:.2f}, '
                     f'minimum coverage of the 2: {min_coverage:.2f}')
        domtab_df = domtab_df[(domtab_df['seq_hmm_len_ratio'] >= min_coverage) &
                              (domtab_df['seq_hmm_len_ratio'] <= 1 / min_coverage) &
                              (domtab_df['hmm_cover'] >= hmm_coverage) &
                              (domtab_df['seq_cover'] >= seq_coverage)].sort_values(
                                  ['seq_name', 'eval'])

        domtab_df = domtab_df.groupby('seq_name', as_index=False).nth(0)
        if output_file is not None:
            SysUtil.mkdir_p(os.path.dirname(output_file))
            domtab_df.to_csv(output_file, index=False)
            logging.info(f'Best matches pared from {domtbl_file} are saved to {output_file}')
        return domtab_df

    @classmethod
    def read_domtab_file_to_df(cls, domtbl_file):
        """
        Parse `domtbl_file`, which is the domtable (--domtblout) output file from
        hmmsearch to data frame, translate the column names to better names.

        Parameters
        ----------
        domtbl_file : str
            path to the domtblout file

        Returns
        -------
        pd.DataFrame
            parsed data frame, with new names for columns, and 3 extra columns
            used for later filtering.
        """
        assert os.path.isfile(domtbl_file), domtbl_file
        col_names = ['seq_name', 'accession_seq', 'seq_len', 'hmm_name', 'accession_hmm', 'hmm_len', 'eval',
                     'full_seq_score', 'full_seq_bias',
                     'domain_idx', 'domain_tot', 'c_eval', 'i_eval', 'domain_score', 'domain_bias', 'hmm_from',
                     'hmm_to', 'seq_from', 'seq_to', 'env_from', 'env_to', 'accuracy']
        # 30 columns in total
        all_col_names = col_names + [f'c_{i}' for i in range(30 - len(col_names))]
        # Skip the first 3 and last 10 rows, a better approach is to remove lines starting with #
        # We cannot use comment='#' option because all the lines contain # character
        # Probably could use fewer columns as well
        domtab_df = pd.read_csv(domtbl_file, sep='\s+', header=None,
                                names=all_col_names,
                                usecols=[col for col in col_names if not col.startswith('accession')],
                                skiprows=3, skipfooter=10, engine='python').astype(
            {'seq_name': str, 'hmm_name': str})
        # calculate the ratios for convenience of filtering later
        domtab_df['seq_hmm_len_ratio'] = domtab_df['seq_len'].values / np.maximum(
            domtab_df['hmm_len'].values, 1e-8)
        domtab_df['hmm_cover'] = (domtab_df['hmm_to'].values -
                                  domtab_df['hmm_from'].values + 1) / domtab_df['hmm_len'].values
        domtab_df['seq_cover'] = (domtab_df['seq_to'].values -
                                  domtab_df['seq_from'].values + 1) / domtab_df['seq_len'].values
        return domtab_df


class ProdigalOutputParser(object):
    '''
    Collection of class methods that parse Prodigal output
    '''
    @classmethod
    def is_header_line(cls, line):
        if not isinstance(line, str):
            return False
        if not (line.startswith('DEFINITION') or line.startswith('# Sequence Data:')):
            return False
        if 'seqnum' not in line:
            return False
        return True

    @classmethod
    def header_to_dict(cls, header_line):
        assert cls.is_header_line(header_line), f'{header_line} is not prodigal header'
        seqnum_start_pos = header_line.find('seqnum')
        return dict(key_val.split('=') for key_val in header_line.strip('\n')[seqnum_start_pos:].split(';'))

    @classmethod
    def cds_note_line_to_dict(cls, line):
        assert isinstance(line, str), line
        line = line.strip().strip('\n')
        assert line.startswith('/note="'), line
        assert line.endswith(';"'), line
        note_content_start_pos = line.find('"')
        return dict(key_val.split('=') for key_val in line[note_content_start_pos + 1:-2].split(';'))

    @classmethod
    def get_seq_type(cls, header_line):
        '''
        From prodigal output header line, infer the sequence type.

        Parameters
        ----------
        header_line : str

        Returns
        -------
        seq_type : {'PLASMID', 'CHRSM', 'CONTIG'}
        plasmid, or complete genome or other (mitochondrial chromosome, chloroplast chromosome)
        '''
        assert isinstance(header_line, str), header_line
        assert '\n' not in header_line, header_line
        if "plasmid" in header_line:
            return "PLASMID"
        if "genome" in header_line or "chromosome" in header_line:
            return "CHRSM"
        return "CONTIG"

    @classmethod
    def is_seq_circular(cls, header_line):
        '''
        Returns whether the sequence is circular/complete from prodigal output headerline

        Returns
        -------
        bool
        '''
        return 'complete' in header_line

    @classmethod
    def get_contig_info(cls, prod_outfile):
        '''
        From prodigal's output (gbk) file, find the information of the contigs
        from which the proteins are called.

        Parameters
        ----------
        prod_outfile: str
            file path to the prodigal output.
            Sample content of such a file
            inGenBank format:
            DEFINITION  seqnum=1;seqlen=4367;...
            seqhdr="gi|226315872|ref|NC_006969.2| Rhodococcus opacus B4 plasmid pKNR01, complete sequence";...
            version=Prodigal.v2.60;run_type=Single;model="Ab initio";gc_cont=67.62;transl_table=11;uses_sd=1
            # in gff format:
            # Sequence Data: seqnum=1;seqlen=5115410;seqhdr="Contig1"


        Returns
        -------
        contigs_info : dict
            key: contig's sequence number in string
            val: the contig's information, including its id, sequence type, circularity,
                 and CDS info for all its proteins
        '''
        contigs_info = dict()
        logging.info(f'Reading prodigal output file {prod_outfile}')
        with open(prod_outfile, 'r') as prod_out:
            for line in prod_out:
                if cls.is_header_line(line):
                    logging.debug(f'New header line {line}, start a new contig, ID_CDS_dict reset to empty')
                    ID_CDS_dict = dict()
                    line = line.strip("\n")
                    seq_key_to_val = cls.header_to_dict(line)
                    contigs_info[seq_key_to_val['seqnum']] = {
                        'seqID': seq_key_to_val['seqhdr'].split()[0],
                        'seqType': cls.get_seq_type(line),
                        'circular': 'Y' if cls.is_seq_circular(line) else 'N',
                        'proteins': ID_CDS_dict}

                elif line.strip().startswith("CDS"):
                    CDS_string = line.strip(' \n').split()[-1]
                    note_line = prod_out.readline()
                    note_dict = cls.cds_note_line_to_dict(note_line.strip('\n'))
                    contig_num, protein_num = note_dict['ID'].split('_')
                    contigs_info[contig_num]['proteins'][protein_num] = CDS_string

                else:
                    continue

        return contigs_info

    @classmethod
    def write_genetic_element(cls, genetic_element_file, contig_info, annot_format="pf"):
        '''
        Write genetic-elements.dat file with the presence of prodigal's output file,
        which provides more information about the contigs/sequences in the genome

        Input:
        1.  genetic_elem - file object to write the content
        2.  contig_info
        3.  annot_format - annotation format, either .pf or .gbk

        Output:
        1.  write content to the genetic-elements.dat
        '''

        SysUtil.mkdir_p(os.path.join(genetic_element_file))
        # counter for each type of sequence
        plasmid_ct = 0
        chrsm_ct = 0
        contig_ct = 0
        contig_to_annot_file = dict()

        with open(genetic_element_file, 'w') as genetic_elem:
            for i in contig_info:
                seq_type = contig_info[i]['seqType']
                if seq_type == "CHRSM":
                    chrsm_ct += 1
                    element_id = seq_type + "-" + str(chrsm_ct)
                elif seq_type == "PLASMID":
                    plasmid_ct += 1
                    element_id = seq_type + "-" + str(plasmid_ct)
                elif seq_type == 'CONTIG':
                    contig_ct += 1
                    element_id = seq_type + "-" + str(contig_ct)
                else:
                    raise ValueError(f'seq_type {seq_type} must be CHRSM, PLASMID, or CONTIG')

                genetic_elem.write(f'ID\t{element_id}\n')
                genetic_elem.write(f'TYPE\t:{seq_type}\n')
                genetic_elem.write(f'CIRCULAR?\t{contig_info[i]["circular"]}\n')

                annot_file = f'{element_id}.{annot_format}'
                contig_to_annot_file[i] = annot_file

                genetic_elem.write(f'ANNOT-FILE\t{annot_file}\n')
                genetic_elem.write(f'SEQ-FILE\t \n')
                genetic_elem.write("//\n")
        logging.info(f'genetic_element info is written to {genetic_element_file}')


class OutputReader(object):
    """
    Class methods to help read output from various programs
    """
    @classmethod
    def read_tRNAscan_output(cls, tRNAscan_output_file):
        """
        Parse tRNAscan's output file *.o from -o option into data frame.

        sample .tRNAscan.o content is below:

        Sequence                tRNA    Bounds  tRNA    Anti    Intron Bounds   Inf
        Name            tRNA #  Begin   End     Type    Codon   Begin   End     Score   Note
        --------        ------  -----   ------  ----    -----   -----   ----    ------  ------
        Contig1         1       1657    1733    Ile     GAT     0       0       73.8
        Contig1         2       1797    1872    Ala     TGC     0       0       81.6

        """
        assert os.path.isfile(tRNAscan_output_file), tRNAscan_output_file
        assert tRNAscan_output_file.endswith('.o'), tRNAscan_output_file
        df = pd.read_csv(tRNAscan_output_file, sep='\t', skiprows=3, header=None,
                         keep_default_na=False, index_col=False,
                         names=['seq_name', 'tRNA_idx', 'tRNA_begin', 'tRNA_end',
                                'tRNA_type', 'anti_codon',
                                'intron_begin', 'intron_end', 'inf_score', 'Note'])
        return df

    @classmethod
    def read_RNAmmer_gff(cls, RNAmmer_gff_file):
        """
        Parse RNAmmer's output file *.gff from -gff option into data frame.

        sample .gff file content is below, the lines not starting with # are separated with \t

        ##gff-version2
        ##source-version RNAmmer-1.2
        ##date 2020-03-30
        ##Type DNA
        # seqname           source                      feature     start      end   score   +/-  frame  attribute
        # ---------------------------------------------------------------------------------------------------------
        Contig1 RNAmmer-1.2     rRNA    2149    5078    3285.0  +       .       23s_rRNA
        Contig1 RNAmmer-1.2     rRNA    7490    10419   3283.7  +       .       23s_rRNA
        Contig1 RNAmmer-1.2     rRNA    5181    5296    77.0    +       .       5s_rRNA
        Contig1 RNAmmer-1.2     rRNA    10513   10628   72.5    +       .       5s_rRNA
        Contig1 RNAmmer-1.2     rRNA    8       1548    1853.3  +       .       16s_rRNA
        Contig1 RNAmmer-1.2     rRNA    5629    7175    1931.7  +       .       16s_rRNA
        # ---------------------------------------------------------------------------------------------------------

        """
        assert os.path.isfile(RNAmmer_gff_file), RNAmmer_gff_file
        assert RNAmmer_gff_file.endswith('.gff'), RNAmmer_gff_file
        df = pd.read_csv(RNAmmer_gff_file, sep='\t', comment='#', header=None,
                         keep_default_na=False, index_col=False,
                         names=['seq_name', 'source', 'feature', 'start', 'end', 'score', 'pm', 'frame', 'attribute'])
        return df
