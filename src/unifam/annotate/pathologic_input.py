'''
Prepare input files for Pathologic input
'''
from unifam.base.util import SysUtil
import logging
import pandas as pd


class PathoLogicInput(object):

    @classmethod
    def prot_cluster_annot_dict_to_patho_dict(cls, cluster_annot_dict):
        '''
        Returns a dictionary with annotations that will be used in pathway tools.
        Some fields are not from annotation, such as STARTBASE and ENDBASE.

        Fields in Pathologic (.pf) file for a gene:
        ABUNDANCE optional
        CODING-SEGMENT optional
        DBLINK optional <DB:Accession> (Use UNIPROT or SP for UniProt accession number)
        EC recommended
        STARTBASE, ENDBASE optional
        FUNCTION required, use "ORF" if unknown
        FUNCTION-SYNONYM optional
        GENE-COMMENT optional
        GO recommended (will be written in DBLINK field)
        ID highly recommended - unique identifier for the gene
        METACYC optional - MetaCyc reaction ID assigned for the protein if known
        NAME required - common name of the gene
        PRODUCT-TYPE required - P, PSEUDO, TRNA, RRNA, MISC-RNA
        SYNONYM

        Example
        -------
        ;;; The PF file format.
        ;;; This is a comment in front of an imaginary example record.
        ;;; It starts with a semicolon. Each record starts with an "ID"
        ;;; line, and is terminated by a "//" line.
        ID b1262
        NAME trpC
        STARTBASE 1317812
        ENDBASE 1316451
        DBLINK SP:P00909
        PRODUCT-TYPE P
        SYNONYM foo
        SYNONYM foo2
        GENE-COMMENT f453; 99 pct identical to TRPC_ECOLI SW:P00909
        ;;; The following shows how information about multiple functions of
        ;;; a protein is supplied:
        FUNCTION N-(5-phosphoribosyl) anthranilate isomerase
        EC 5.3.1.24
        FUNCTION-SYNONYM phosphoribosyl anthranilate isomerase
        FUNCTION-COMMENT Amino acid biosynthesis: Tryptophan (3rd step)
        FUNCTION indole-3-glycerolphosphate synthetase
        EC 4.1.1.48
        FUNCTION-COMMENT	Amino acid biosynthesis: Tryptophan (4th step)
        DBLINK GO:0000250
        //
        '''
        patho_key_set = {'GO', 'ec', 'full_name', 'short_name', 'synonym', 'ORF', 'OLN'}
        assert patho_key_set < set(cluster_annot_dict), \
                    f'cluster_annot_dict missing keys {patho_key_set - set(cluster_annot_dict)}'
        patho_dict = dict()
        patho_dict['DBLINK'] = cluster_annot_dict['GO']

        patho_dict['EC'] = cluster_annot_dict['ec']
        patho_dict['FUNCTION'] = cluster_annot_dict['full_name'] + cluster_annot_dict['short_name']
        if not len(patho_dict['FUNCTION']):
            patho_dict['FUNCTION'] = ['ORF']
        patho_dict['FUNCTION-SYNONYM'] = cluster_annot_dict['synonym']
        patho_dict['FUNCTION-COMMENT'] = cluster_annot_dict['ORF'] + cluster_annot_dict['OLN']
        patho_dict['PRODUCT-TYPE'] = ['P']

        return patho_dict

    @classmethod
    def get_prot_patho_dict_list(cls, domtab_df, contig_to_info, cluster_annot_dict):
        """
        From hmmsearch result `domtab_df` (protein --> hmm),
        prodigal output gbk result `contig_to_info` (protein --> start, end bases)
        `cluster_annot_dict`, generate protein annotations in pathologic dict format
        """
        assert isinstance(domtab_df, pd.DataFrame), type(domtab_df)

    @classmethod
    def rrna_df_to_patho_dict(cls, rrna_series, seq_name_to_id):
        """
        Convert RNAmmer output parse result df to pathologic key,val pair dict.
        This dict will then translate to a record in the annotation .pf file

                                        seq_name       source feature  start  end  score pm frame attribute
        0  NODE_11820_length_3376_cov_428.234869  RNAmmer-1.2    rRNA     26  996  154.1  +     .  16s_rRNA

        """
        assert isinstance(rrna_series, pd.Series), type(rrna_series)
        assert isinstance(seq_name_to_id, dict), type(seq_name_to_id)
        patho_dict = dict()
        seq_name = rrna_series['seq_name']
        patho_dict['NAME'] = seq_name_to_id[seq_name]
        patho_dict['ID'] = f'{seq_name}_rRNA_{rrna_series["attribute"]}'
        strand_pm = rrna_series['pm']
        start_base = rrna_series['start'] if strand_pm == '+' else rrna_series['end']
        end_base = rrna_series['end'] if strand_pm == '+' else rrna_series['start']
        patho_dict['STARTBASE'] = str(start_base)
        patho_dict['ENDBASE'] = str(end_base)
        patho_dict['PRODUCT-TYPE'] = ['RRNA']
        return patho_dict

    @classmethod
    def trna_df_to_patho_dict(cls, trna_series, seq_name_to_id):
        """
        Convert tRNAscan output parse result df to pathologic key,val pair dict.
        This dict will then translate to a record in the annotation .pf file

                                 seq_name  tRNA_idx  tRNA_begin  tRNA_end tRNA_type anti_codon  intron_begin  intron_end
0     NODE_1_length_526561_cov_50.907785          1      120500    120584       Ser        TGA             0           0
1     NODE_1_length_526561_cov_50.907785          2      202054    202127       His        GTG             0           0
2     NODE_1_length_526561_cov_50.907785          3      245494    245581       Ser        GGA             0           0
3     NODE_1_length_526561_cov_50.907785          4      523248    523332       Ser        TGA             0           0
        """
        assert isinstance(trna_series, pd.Series), type(trna_series)
        assert isinstance(seq_name_to_id, dict), type(seq_name_to_id)
        patho_dict = dict()
        seq_name = trna_series['seq_name'].strip()
        patho_dict['NAME'] = f"{seq_name_to_id[seq_name]}_{trna_series['tRNA_idx']}"
        patho_dict['ID'] = f"{seq_name}_tRNA_idx_{trna_series['tRNA_idx']}"
        patho_dict['STARTBASE'] = str(trna_series['tRNA_begin'])
        patho_dict['ENDBASE'] = str(trna_series['tRNA_end'])
        patho_dict['PRODUCT-TYPE'] = ['TRNA']
        return patho_dict

    @classmethod
    def pf_str_from_patho_dict(cls, patho_dict):
        """
        Given a dict containing protein annotation information,
        convert it to string in pathologic input format.

        The keys in `patho_dict` are the annotation fields in pathologic input.
        """
        assert isinstance(patho_dict, dict)
        pf_str_list = []
        for key, val in patho_dict.items():
            if isinstance(val, str):
                pf_str_list.append(f'{key}\t{val}')
            elif isinstance(val, list):
                for val_element in val:
                    pf_str_list.append(f'{key}\t{val_element}')
            else:
                raise TypeError(f'{type(val)} is not str or list')
        pf_str_list.append('//')
        return '\n'.join(pf_str_list)

    @classmethod
    def write_organism_params(cls, file_path, org_param_dict):
        SysUtil.mkdir_p(file_path)
        sep = '\t'
        with open(file_path, 'w') as org_param_f:
            for key, val in org_param_dict.items():
                org_param_f.write(f'{key}{sep}{val}\n')
        logging.info(f'organism params data written to {file_path}')

    @classmethod
    def get_organism_params_dict(cls, organism, db_id, domain=None,
                                codon_table='11', ncbi_tax_id=None):
        '''
        organism: config
        db_id: 'U' + prefix.upper()
        codon_table: prodigal gffFile
        ncbi_tax_id: config
        '''
        DOMAIN_TO_TAXID = {'bac': 2, 'arc': 2157, 'euk': 2759, 'unknown': 131567}
        if domain is None:
            domain = 'unknown'
        assert domain in DOMAIN_TO_TAXID, DOMAIN_TO_TAXID

        # 2 to 10 alphanumeric characters, no intervening spaces, must start with a letter
        org_param_dict = {'ID': db_id,
                          'STORAGE': 'FILE',
                          # full species name, e.g., Escherichia coli
                          'NAME': organism,
                          'PRIVATE?': 'NIL',
                          'RANK': '|strain|',
                          'CODON-TABLE': codon_table,
                          'DBNAME': f'{db_id}Cyc',
                          'CREATE?': 't',
                          'DOMAIN': f'TAX-{DOMAIN_TO_TAXID[domain]}',
                         }
        if ncbi_tax_id is not None:
            org_param_dict['NCBI-TAXON-ID'] = ncbi_tax_id

        return org_param_dict
