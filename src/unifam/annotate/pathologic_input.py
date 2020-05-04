'''
Prepare input files for Pathologic input
'''
from unifam.base.util import SysUtil
import logging


class PathoLogicInput(object):

    @classmethod
    def pf_str_from_patho_dict(cls, patho_dict):
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
