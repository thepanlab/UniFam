#!/usr/bin/python

'''Parse the .uc output from usearch, which includes information about every cluster
Output format: seq_name [tab] cluster_number
'''
import sys
from subprocess import Popen, PIPE, check_call, STDOUT
from datetime import datetime
## =================================================================
## argument parser
import argparse


## =================================================================
def parse_uc(ucfile,group_sg,count=None,write_count=False):
    ''' Parse .uc file output from usearch clust, to get simpler summary of the groups.
        Input:  .uc file output from usearch clust
        Output: group_sg - tab delimited file that has 3 columns: sequenceID, group number, seed length
                count - file with count of sequences in each group
    '''
    group_sg.write("{}\t{}\t{}\n".format("sequenceID", "groupID", "seed_length"))
    with open(ucfile,'r') as f:
        for line in f:
            line = line.strip("\n")
            line = line.split()
            if line[0] != 'C':
                ## sequenceID \t cluster_number \t length of seed sequence
                ## if the sequence is not seed, write * for the length
                group_sg.write('{0:}\t{1:}\t{2:}\n'.format(line[-2],line[1], line[2] + '*' if line[-1]=='*' else line[2]))
            else: 
                if write_count:
                    if count is None:
                        count = sys.stdout
                    count.write(line[1] + "\t" + line[2]+'\n') ## number of proteins in each group: groupID count

def sg_to_gs(sg_file,gs_out):
    ''' from group file of format [seq group] to format [group : seq : count : seed seq length]
    '''
    group_dict = dict()
    with open(sg_file, 'r') as sg:
        sg.readline()
        for line in sg:
            seqID, groupID, length = line.strip('\n').split('\t')
            groupID = int(groupID)
            if groupID in group_dict:
                group_dict[groupID].append(seqID)
            else:
                group_dict[groupID] = [seqID]
    for group in sorted(group_dict):
        gs_out.write('{}\t{}\n'.format(group, ' '.join(group_dict[group])))

parser = argparse.ArgumentParser(description="Parse .uc file to find the groups and statistics of the groups",
    prog = 'parse_uc', #program name
    prefix_chars='-', # prefix for options
    fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
    conflict_handler='resolve', # for handling conflict options
    add_help=True, # include help in the options
    formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
    )

## input group file
parser.add_argument("ucfile",help=".uc file to parse")
## optional output (count sequences in each group)
parser.add_argument("--count",action='store_true',help="save the count file?",dest="write_count")
parser.add_argument("-c", "--countfile",help="file name for group counts",dest='countfile',default=sys.stdout, type=argparse.FileType('w'))
## output file name
parser.add_argument("-g","--group",help="group file with seq -- group format",dest='groupfile')

## ===========
## Main function
## ===========
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = parser.parse_args(argv)
    if args.groupfile is not None:
        group_sg_out = open(args.groupfile+'_sg','w')
        group_gs_out = open(args.groupfile+'_gs','w')
    else:
        group_sg_out = sys.stdout
        group_gs_out = sys.stdout
    start_time = datetime.now()
    sys.stderr.write("\n===============================================================================\n")
    sys.stderr.write("parse_uc starts running\n")
    ## parse the uc file
    parse_uc(args.ucfile,group_sg_out,args.countfile,args.write_count)
    if args.groupfile is not None:
        group_sg_out.close()
        ## get the group_seq file from seq_group file
        sg_to_gs(args.groupfile+'_sg',group_gs_out)
        group_gs_out.close()
    
    #sys.stderr.write("sort the group files\n")
    ### order seq_group file by seqID
    #p = Popen("sort " + args.groupfile + "_sg -o " + args.groupfile + "_sg",shell=True,stdout=PIPE,stderr=STDOUT)
    #p.wait()
    ### order group_seq file by group size
    #p = Popen("sort -t: -nk2,2 -nk4,4 " + args.groupfile + "_gs -o " + args.groupfile + "_gs",shell=True,stdout=PIPE,stderr=STDOUT)
    #p.wait()
    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write("Total elapsed time is %s [seconds]\n" %duration)
    sys.stderr.write("\n===============================================================================\n")
## ======================    
## call from command line
## ======================    
if __name__ == "__main__":
    sys.exit(main())
