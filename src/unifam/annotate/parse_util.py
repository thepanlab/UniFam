'''
Utility classes for parsing output
'''
import logging


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
        From prodigal's output file, find the information of the contigs from which the proteins are called.
        Supposedly it should work for both gff format (default for prodigal3) and gbk format

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


