'''
Prepare input files for Pathologic input
'''
from unifam.base.util import SysUtil
import logging
import pandas as pd
import os


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
    def write_genetic_element(cls, genetic_element_file,
                              contig_to_info, domtab_df, cluster_annot_dict,
                              rrna_df, trna_df):
        '''
        Write genetic-elements.dat file with the presence of prodigal's output file,
        which provides more information about the contigs/sequences in the genome

        Input:
        1.  genetic_elem - file object to write the content
        2.  contig_to_info

        Output:
        1.  write content to the genetic-elements.dat
        '''
        patho_dir = os.path.dirname(genetic_element_file)
        SysUtil.mkdir_p(patho_dir)
        # counter for each type of sequence
        plasmid_ct = 0
        chrsm_ct = 0
        contig_ct = 0
        contig_to_annot_file = dict()

        with open(genetic_element_file, 'w') as genetic_elem:
            for i, contig_info in contig_to_info.items():
                seq_type = contig_to_info[i]['seqType']
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

                seq_id = contig_to_info[i]['seqID']
                genetic_elem.write(f'ID\t{element_id}\n')
                genetic_elem.write(f'NAME\t{seq_id}\n')
                genetic_elem.write(f'TYPE\t:{seq_type}\n')
                genetic_elem.write(f'CIRCULAR?\t{contig_info["circular"]}\n')

                annot_file = f'{element_id}.pf'
                contig_to_annot_file[i] = annot_file
                annot_file_path = os.path.join(patho_dir, annot_file)

                cls._write_rna_annot_file(annot_file_path, rrna_df, seq_id, 'rrna', 'w')
                cls._write_rna_annot_file(annot_file_path, trna_df, seq_id, 'trna', 'w')
                cls._write_protein_annot_file(annot_file_path, i, contig_info,
                                              domtab_df, cluster_annot_dict, 'a')

                genetic_elem.write(f'ANNOT-FILE\t{annot_file}\n')
                genetic_elem.write(f'SEQ-FILE\t \n')
                genetic_elem.write("//\n")
                print(f'{seq_id} writtein to {annot_file}')
        logging.info(f'genetic_element info is written to {genetic_element_file}')

    @classmethod
    def _write_protein_annot_file(cls, annot_file, contig_num, contig_info,
                                  domtab_df, cluster_annot_dict, file_mode='a'):
        """
        Write protein annotations for a contig, with the following informatin:
        * `contig_num` contig_info is an item pair from contig_to_info parsed from gbk file
        * `domtab_df` provides the mapping from protein sequence name to best matching hmm name
        * `cluster_annot_dict` provides the mapping from hmm name to the corresponding cluster's annotation
        """
        assert isinstance(domtab_df, pd.DataFrame), type(domtab_df)
        seq_id = contig_info['seqID']
        with open(annot_file, file_mode) as annot_pf:
            for protein_num, cds in contig_info['proteins'].items():
                protein_seq_name = f'{seq_id}_{protein_num}'
                start_pos, end_pos = cds
                protein_subdf = domtab_df[domtab_df['seq_name'] == protein_seq_name]
                if protein_subdf.empty:
                    continue
                assert len(protein_subdf) == 1, (protein_seq_name, protein_subdf)
                protein_hmm = protein_subdf['hmm_name'].values[0]
                protein_patho_dict = cls.prot_cluster_annot_dict_to_patho_dict(
                    cluster_annot_dict[protein_hmm])
                protein_patho_dict['STARTPOS'] = str(start_pos)
                protein_patho_dict['ENDPOS'] = str(end_pos)
                protein_patho_dict['NAME'] = protein_seq_name
                protein_patho_dict['ID'] = f'{contig_num}_protein_{protein_num}'
                annot_pf.write(cls.pf_str_from_patho_dict(protein_patho_dict))

    @classmethod
    def _write_rna_annot_file(cls, annot_file, rna_df, seq_name, rna_type, file_mode='a'):
        assert rna_type in ('trna', 'rrna'), rna_type
        seq_sub_df = rna_df[rna_df['seq_name'] == seq_name]
        if not seq_sub_df.empty:
            with open(annot_file, file_mode) as annot_pf:
                for _, r in seq_sub_df.iterrows():
                    patho_dict = PathoLogicInput.rrna_df_to_patho_dict(r) if rna_type == 'rrna' \
                            else PathoLogicInput.trna_df_to_patho_dict(r)
                    annot_pf.write(PathoLogicInput.pf_str_from_patho_dict(patho_dict))

    @classmethod
    def rrna_df_to_patho_dict(cls, rrna_series):
        """
        Convert RNAmmer output parse result df to pathologic key,val pair dict.
        This dict will then translate to a record in the annotation .pf file

                                        seq_name       source feature  start  end  score pm frame attribute
        0  NODE_11820_length_3376_cov_428.234869  RNAmmer-1.2    rRNA     26  996  154.1  +     .  16s_rRNA

        """
        assert isinstance(rrna_series, pd.Series), type(rrna_series)
        patho_dict = dict()
        seq_name = rrna_series['seq_name']
        rrna_idx = rrna_series['rRNA_idx']
        patho_dict['NAME'] = f'rRNA_{rrna_idx}'
        patho_dict['ID'] = f'{seq_name}_rRNA_{rrna_series["attribute"]}'
        strand_pm = rrna_series['pm']
        start_base = rrna_series['start'] if strand_pm == '+' else rrna_series['end']
        end_base = rrna_series['end'] if strand_pm == '+' else rrna_series['start']
        patho_dict['STARTBASE'] = str(start_base)
        patho_dict['ENDBASE'] = str(end_base)
        patho_dict['PRODUCT-TYPE'] = ['RRNA']
        return patho_dict

    @classmethod
    def trna_df_to_patho_dict(cls, trna_series):
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
        patho_dict = dict()
        seq_name = trna_series['seq_name'].strip()
        patho_dict['NAME'] = f"tRNA_{trna_series['tRNA_idx']}"
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
        pf_str_list.append('//\n')
        return '\n'.join(pf_str_list)

    @classmethod
    def write_organism_params(cls, file_path, org_param_dict):
        SysUtil.mkdir_p(os.path.dirname(file_path))
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
