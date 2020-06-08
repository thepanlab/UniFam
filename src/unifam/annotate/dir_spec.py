'''
directory and file structures.
'''
from configparser import ConfigParser
import os
import re
import uuid


class UniFamDirSpec(object):
    """
    Specification of the directory structure for unifam pipeline

    Parameters
    ----------
    work_dir : str
        root of the working directory, where all the intermediate and final output is saved
    prefix : str
        the files/directories will have this string as prefix
    """

    def __init__(self, work_dir, prefix=None):
        assert isinstance(work_dir, str), work_dir
        if prefix is None:
            prefix = str(uuid.uuid4())

        assert isinstance(prefix, str), prefix
        assert prefix, 'prefix is empty'
        # replace any non-alphanumeric character or underscore to underscore
        self._prefix = re.sub(r'\W', r'_', prefix)
        self._work_dir = os.path.expandvars(work_dir)

    @classmethod
    def from_config(cls, config):
        assert isinstance(config, ConfigParser), type(config)
        prefix = config.get('UniFam', 'name')
        work_dir = config.get('UniFam', 'work_dir')
        return cls(work_dir, prefix)

    def get_path_prefix(self):
        return os.path.join(self._work_dir, self._prefix)

    def get_protein_annot_file(self, use_faa_format=False):
        """
        Returns the file path under work_dir for protein annotation file.

        Parameters
        ----------
        use_faa_format : bool
            if True, returns file name for faa with annotation at the header lines
            otherwise returns file name for a flat file with annotation for each protein

        Returns
        -------
        file_path : str
        """
        assert isinstance(use_faa_format, bool)
        if use_faa_format:
            return os.path.join(self._work_dir, f'{self._prefix}_annot.faa')
        return os.path.join(self._work_dir, f'{self._prefix}.annot')

    def get_readme_file(self):
        return os.path.join(self._work_dir, 'README')

    def get_prodigal_out_file(self, out_type):
        assert out_type in {'faa', 'gbk', 'out'}
        return os.path.join(self._work_dir, f'{self._prefix}.prod.{out_type}')

    def get_RNAmmer_gff_file(self):
        return os.path.join(self._work_dir, f'{self._prefix}.RNAmmer.gff')

    @classmethod
    def get_tRNAscan_out_file(cls, out_prefix, out_type):
        assert isinstance(out_prefix, str), out_prefix
        assert out_type in {'o', 'structure', 'stat'}, out_type
        return f'{out_prefix}.tRNA.{out_type}'

    def get_tRNAscan_output_file(self):
        return self.get_tRNAscan_out_file(self.get_path_prefix(), 'o')

    def get_hmmsearch_domtbl_file(self):
        return os.path.join(self._work_dir, f'{self._prefix}.hmmsearch.domtbl')

    def get_hmmsearch_out_file(self):
        return os.path.join(self._work_dir, f'{self._prefix}.hmmsearch.out')

    def get_patho_input_dir(self):
        """
        Returns the directory that contains the input files for pathway-tools.
        """
        return os.path.join(self._work_dir, f'{self._prefix}_patho_input')
