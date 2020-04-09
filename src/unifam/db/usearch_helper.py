'''
This module contains utilities to parse usearch results.
'''
import os
import sys


class UsearchHelper(object):

    @classmethod
    def parse_cluster_file(cls, uc_file):
        '''
        Parse usearch cluster (uc) format file.
        See reference on this page: http://www.drive5.com/usearch/manual/opt_uc.html
        '''
        assert os.path.isfile(uc_file), uc_file

    @classmethod
    def parse_uc(cls, ucfile, group_sg_file, count_file=None, write_count=False):
        """ Parse .uc file output from 'usearch clust' to get simpler summary of the groups.
        Parameters
        -----
        ucfile : str
            .uc file output from 'usearch clust'
        group_sg_file : str
            tab delimited group file that has 3 columns: sequenceID, group number, seed length
        count_file : str
            file to save the count of sequences in each group
        """
        group_sg_file.write("{}\t{}\t{}\n".format(
            "sequenceID", "groupID", "seed_length"))
        with open(ucfile, 'r') as f:
            for line in f:
                line = line.strip("\n")
                line = line.split()
                if line[0] != 'C':
                    # sequenceID \t cluster_number \t length of seed sequence
                    # if the sequence is not seed, write * for the length
                    group_sg_file.write('{0:}\t{1:}\t{2:}\n'.format(
                        line[-2], line[1], line[2] + '*' if line[-1] == '*' else line[2]))
                else:
                    if write_count:
                        if count_file is None:
                            count_file = sys.stdout
                        # number of proteins in each group: groupID count
                        count_file.write(line[1] + "\t" + line[2] + '\n')
