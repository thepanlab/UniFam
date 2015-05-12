#!/usr/bin/python

'''
annotation

Assign annotations to each cluster of protein sequences.
EC name, Recommended full name, Go Term, OLN, ORF;
To do: Keywords (excludes the technical terms and such)

3 ways to do the annotation: 1. intersection; 2. union; 3. consensus

Created by JJ Chai on 11/21/2013
Last modified 12/05/2013
Copyright (c) 2013 JJ Chai (ORNL). All rights reserved.
 
'''
# Import Python modules
import sys, warnings, os
import time
import glob
from datetime import datetime
import re
#from subprocess import Popen, PIPE, check_call, STDOUT

# Import local modules
sys.path.append("/chongle/jj/01_PfClust/scripts")
sys.path.append("/Users/cjg/Work/01_PfClust/scripts")
import jj_utils

## =================================================================
## argument parser
## =================================================================
import argparse

## Version
version_str = "0.0.1"
'''first version of annotation assignment
'''

parser = argparse.ArgumentParser(description="Assign annotations for all the clusters",
    prog = 'annotation', #program name
    prefix_chars='-', # prefix for options
    fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
    conflict_handler='resolve', # for handling conflict options
    add_help=True, # include help in the options
    formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
    )
## version control
parser.add_argument("--version", action="version",version='%(prog)s {}'.format(version_str))

## --verbose mode, default: false (quiet mode)
parser.add_argument("-v", "--verbose", action="store_true",help="verbose mode, more output")

## method to retrieve annotation for group
parser.add_argument("--int", action="store_true",help="annotation by intersection")
parser.add_argument("--con", action="store_true",help="annotation by consensus")

## input sequence annotation file and the group file
parser.add_argument("-t","--table",help="table with annotations for all the protein sequences",dest='tabfile',required=True)
parser.add_argument("-g","--group",help="group file with all the sequences and their group membership",dest='groupfile',required=True)
## output annotation file
parser.add_argument("-o", "--out",nargs='?',help="output file of the cluster annotations",const='anno.txt',default='anno.txt',dest='outputfile')
## =================================================================
## annotation function
## =================================================================

def annotation(tabfile,groupfile,outputfile,method):
    
    '''Assign annotations to each cluster of protein sequences.
        EC name, Recommended full name, Go Term, OLN, ORF;
    '''
    ## dictionary with groupID as key and seqIDs as the mapped value for the corresponding key
    ## groupID -> list(seqIDs), find for each cluster all the sequences in it
    group_dict = dict()
    with open(groupfile, 'r') as group:
        for line in group: # format: sp|A0A183|LCE6A_HUMAN   243797  2.40375
            group_info = line.strip("\n").split()
            seqID = group_info[0]
            if seqID[0] == 's':
                seqID = seqID.split('|')[1] # only the accession number (primary ID) is extracted, e.g. A0A183
                groupID = int(group_info[1]) # groupID as integer, so that it can be ordered (by number, i.e. length)
                # if groupID is not already a key in the dictionary, create
                if( groupID not in group_dict):
                    group_dict[groupID] = [seqID];
                # otherwise, append the sequenceID to the list of sequences for this group
                else:
                    group_dict[groupID].append(seqID);
            else:
                continue

    ## echo the number of groups to annotate
    print "number of groups is ", len(group_dict)

    ## sequence dictionary
    seq_dict = dict()
#    anno_keys = ['org', 'species','ECname', 'gene_name', 'OLN', 'ORF', 'GO', 'KW', 'full_name']
    ## species names are not included, just organism name
    anno_keys = ['org','ECname', 'gene_name', 'OLN', 'ORF', 'full_name', 'GO', 'KW' ]
    with open(tabfile,'r') as table:
        next(table) # skip the first line (header)
        for line in table:
            line_list = line.strip("\n").split("\t") # split into list
            AC = line_list[0].strip() # accession number (primary id)
            organism = line_list[1].strip() # Organism
#            species = line_list[2].strip() # Species
            ECname = line_list[3].strip() # Recommended EC number
            ECname_set = set(ECname.split(":")) # convert to list

            gene_name = line_list[4].strip().upper()
            gene_name_set = set(gene_name.split(":")) # list of gene names (1 name for each protein/?)

            OLN = line_list[5].strip()
            OLN_set = set(OLN.split(":")) # Ordered Locus Names

            ORF = line_list[6].strip() # temporary ORF names
            ORF_set = set(ORF.split(":"))

            go_terms = line_list[7].strip() # GO terms
            go_terms_set = set(go_terms.split(":"))

            keywords = line_list[8].strip() # Keywords
            keywords_set = set(keywords.split(":"))

            full_name = line_list[9].strip()# Recommended full names, or product name, converted to lower case
            full_name_set = set(full_name.split(":"))

            ## create set for each field of annnotation, later will take union, or intersection, or consensus
            seq_dict[AC] = {'org' : set([organism]), 'ECname' : set(ECname_set), 'gene_name' : gene_name_set,'OLN' : OLN_set, 'ORF' : ORF_set, 'GO': go_terms_set, 'KW': keywords_set, 'full_name': full_name_set}
    ## echo the number of sequences in the tab file
    print "number of sequence annotations read is ", len(seq_dict)

    ## write output: annotation for each group
    with open(outputfile,'w') as output:
        output.write( "groupID" + "\t" + "size"  + "\t" + "\t".join(anno_keys) + "\n")
        ## loop through all keys
        for key in group_dict:

            full_name_count = 1 # conflict of full name
            group_size = len(group_dict[key]) # number of sequences in the group
            output.write(str(key) + "\t" + str(group_size) + "\t") # groupID group_size
            
            ## loop through all annotation fields
            for anno_key in anno_keys:
                # initialize the set of annotation for this annotation field
                # from the first sequence
                anno_set = seq_dict[group_dict[key][0]][anno_key]
                anno_name = " ".join(anno_set); # as a string
#                anno_name = re.sub('[^0-9a-z]',' ', anno_name) # replace non alphanumeric symbols by blanks

                ## loop through all the rest of sequences in the group
                for i in xrange(1, group_size):
                    new_anno_set = seq_dict[group_dict[key][i]][anno_key] # get annotation for the next sequence in the group

                    ## if next sequence has a different annotation, add to the list
                    ## union method to annotate
                    if len(new_anno_set.difference(anno_set)) >= 1:
#                        full_name_count += 1
                        new_anno_name = " ".join(new_anno_set) # new annotation phrase
#                        new_anno_name = re.sub('[^0-9a-z]',' ',new_anno_name)
#                        anno_name = jj_utils.lcs(anno_name.split(), new_anno_name.split()) # find the lcs of the annotations
                        anno_set.update(new_anno_set) # update the annotation set
            
                ## delete NA in the set if it's not the only available annotation
                if len(anno_set) > 1:
                    anno_set.discard("NA")
                # write this annotation field to file
                output.write(":".join(anno_set) + "\t")

            output.write("\n")


## =================================================================
## main function
## =================================================================
def main(argv=None):

    if argv is None:
        args = parser.parse_args()

    if args.int:
        method = "intersect"
    elif args.con:
        method = "consensus"
    else:
        method = "union"
    ## print some information
    if args.verbose:
        print 'running verbosely'
        print 'input table file is: {}'.format(args.tabfile)
        print 'input group file is: {}'.format(args.groupfile)
        print 'output file is: {}'.format(args.outputfile)
        print 'annotation method is: {}'.format(method)
    else:
        print 'running quietly\n'
    ## display work start, and time record
    start_time = datetime.now()
    sys.stderr.write("\n===============================================================================\n")
    sys.stderr.write("Assigning annotations for each group\n")
    ## run program and generate the fasta files
    annotation(args.tabfile,args.groupfile,args.outputfile,method)
    ## display work end, and time record, elapsed time
    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write("\nTotal Elapsed Time = [%s] [seconds] \n" % jj_utils.format_time(duration))
    sys.stderr.write("===============================================================================\n")

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
