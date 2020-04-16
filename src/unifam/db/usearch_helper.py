'''
This module contains utilities to parse usearch results.
'''
import os
import pandas as pd


class UsearchHelper(object):
    '''
    TODO: rename the class

    Usearch and annotation related methods.

    Parameters
    ----------
    uc_file : str
    '''

    def __init__(self, uc_file):
        self._uc_df = self.read_usearch_cluster_file_to_df(uc_file)

    def get_sp_cluster_set(self):
        '''
        Returns
        -------
        set
            set of cluster numbers that contain swissprot sequences.
        '''
        return set(self._uc_df[(self._uc_df['record_type'] != 'C') &
                               (self._uc_df['query_label'].str.startswith('sp'))]['cluster_number'].values)

    def get_seqs_in_cluster(self, cluster_number):
        '''
        Parameters
        ----------
        cluster_number : int
            cluster number (index)

        Returns
        -------
        list of str
            list of sequence names in the cluster with given `cluster_number`.
        '''
        cluster_df = self._uc_df[self._uc_df['cluster_number'] == cluster_number]
        seq_list = cluster_df[cluster_df['record_type'] != 'C']['query_label'].values.tolist()
        num_seq_from_C = cluster_df[cluster_df['record_type'] == 'C']['length'].values[0]
        assert len(seq_list) == num_seq_from_C, (f'number of sequences {len(seq_list)} != '
                                                 f'cluster size from C records {num_seq_from_C}')
        return seq_list

    @classmethod
    def read_usearch_cluster_file_to_df(cls, uc_file):
        '''
        Read tab separrated usearch cluster (uc) format file.
        See reference on this page: http://www.drive5.com/usearch/manual/opt_uc.html

        Parameters
        ----------
        uc_file : str
            path to the usearch cluster file

        Returns
        -------
        pd.DataFrame
            uc file pared to data frame, with dtype properly set,
            and the unused columns removed.
        '''
        assert os.path.isfile(uc_file), uc_file

        # record_type H: hit, C: cluster record, S: centroid
        # C records length field is the size of the cluster
        # H,S records length field is the length of the query sequence
        field_names = ['record_type', 'cluster_number', 'length', 'pct_identity',
                       'strand', 'not_used_0', 'not_used_1', 'alignment', 'query_label', 'target_label']
        use_cols = [col_name for col_name in field_names if not col_name.startswith('not_used_')]

        dtype_dict = {col_name: str for col_name in use_cols}
        dtype_dict.update({'cluster_number': int, 'length': int, 'pct_identity': float})

        uc_df = pd.read_csv(uc_file, sep='\t', header=None, names=field_names, usecols=use_cols,
                            keep_default_na=False, na_values={'pct_identity': '*'},
                            dtype=dtype_dict)
        return uc_df
