#!/usr/bin/python

'''
annotation

Assign annotations to each cluster of protein sequences.
EC name, Recommended full name, Go Term, OLN, ORF;

Created by JJ Chai on 11/21/2013
Last modified 05/20/2015
Copyright (c) 2013 JJ Chai (ORNL). All rights reserved.
 
'''
# Import Python modules
import sys, warnings, os
import time
import glob
from datetime import datetime
import re
from Bio import SwissProt

## =================================================================
## argument parser
## =================================================================
import argparse

## Version
version_str = "0.0.2"
''' Second version of annotation assignment.
    Changes:
        Using Bio.SwissProt (read.GN is modified to conform to the new format with {evidence}).
        Add taxonomy lineage information for each group, cleaner code.
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

## input sequence annotation file and the group file
parser.add_argument("-d","--dat",help="SwissProt annotation dat file",dest='dat_file',required=True)
parser.add_argument("-g","--group",help="group file with all the sequences and their group membership",dest='groupfile',required=True)
## output annotation file
parser.add_argument("-o", "--out",nargs='?',help="output file of the cluster annotations",const='annot.txt',default='annot.txt',dest='outputfile')
## =================================================================
## annotation function
## =================================================================

def annotation(sprot_dat_file,groupfile,outputfile):
    ''' Assign annotations to each cluster of protein sequences.
        EC name, Recommended full name, Go Term, OLN, ORF;
    '''
    ## groupID -> list(seqIDs), find for each cluster all the sequences in it
    group_dict, seq_dict = get_group_dict(groupfile)
    seq_dict = read_sprot_dat(sprot_dat_file, seq_dict)
    annot_keys = ['gene_name', 'EC', 'full_name', 'GO', 'organism', 'tax_id', 'lineage', 'OLN', 'ORF', 'KW'] # map primary ID to annotation dictionary

    ## write output: annotation for each group
    with open(outputfile,'w') as output:
        output.write( "groupID" + "\t" + "size"  + "\t" + "\t".join(annot_keys) + "\n")
        ## loop through all groups
        for key in group_dict:
            group_size = len(group_dict[key]) # number of sequences in the group
            output.write("{}\t{}".format(key, str(group_size))) # groupID group_size
            gannot_dict = group_annot(group_dict[key], seq_dict)
            for annot_key in annot_keys:
                my_list = gannot_dict[annot_key]
                if len(my_list) > 0:
                    output.write("\t{}".format(":".join(my_list)))
                else:
                    output.write("\tNA")
            output.write("\n")

def group_annot(seqIDs, seq_dict):
    ''' Get common annotation for a list of sequence IDs
    '''
    gannot_dict = dict()
    myseq_dict = {seqID:seq_dict[seqID] for seqID in seqIDs} # smaller dictionary only for the sequences in this group, to reduce the giant dictionary query
    ID0 = seqIDs[0]
    annot_keys = ['gene_name', 'EC', 'full_name', 'GO', 'organism', 'tax_id', 'lineage', 'OLN', 'ORF', 'KW'] # map primary ID to annotation dictionary
    for key in annot_keys: # loop through all the annotation keys
        value = myseq_dict[ID0][key] # dictionary value for this key at the first sequence, as a template
        if key != "lineage": # lineage uses intersection, others use union
            # for string-value keys
            if type(value) is str:
                gannot_dict[key] = list(set([ myseq_dict[seqID][key] for seqID in seqIDs]))
            # for list-value keys, merge all the lists, and take the set
            elif type(value) is list:
                gannot_dict[key] = list(set(sum([myseq_dict[seqID][key] for seqID in seqIDs],[])))
            # should not have other type of values
            else:
                sys.stderr.write("annotation key {} is neither string nor list\n".format(key))
        else: # lineage should be intersection
            g_lineage = list(set.intersection(*(list(map(set, [myseq_dict[seqID][key] for seqID in seqIDs])))))
            sorted_lineage = value
            lineage_indices = [sorted_lineage.index(x) for x in g_lineage]
            gannot_dict[key] = [x for (y,x) in sorted(zip(lineage_indices, g_lineage))]
    return gannot_dict

def tab_sprot_dat(seq_dict, tab):
    ''' Write annotation for SwissProt proteins into a table.
    '''
    annot_keys = ['gene_name', 'EC', 'full_name', 'GO', 'organism', 'tax_id', 'lineage', 'OLN', 'ORF', 'KW'] # map primary ID to annotation dictionary
    tab.write("seqID" + "\t" + "\t".join(annot_keys) + "\n")
    num_record = 0
    for key in seq_dict:
        num_record += 1
        if num_record % 100000 == 0:
            sys.stderr.write("{} records written\n".format(num_record))
        a = seq_dict[key]
        tab.write(key)
        for annot_key in annot_keys:
            a_list = a[annot_key]
            if type(a_list) is list:
                if len(a_list) > 0:
                    tab.write("\t{}".format(":".join(a_list)))
                else:
                    tab.write("\tNA")
            elif type(a_list) is str:
                tab.write("\t" + a_list)
            tab.write("\n")

def get_group_dict(groupfile):
    ''' Read group file with format (seqID\tgroupID), and return two dictionaries:
        group_dict: groupID -> list of seqIDs in this group
        seq_dict: seqID -> {} (empty dictionary to be filled in when reading the sprot_dat file
    '''
    group_dict = dict()
    seq_dict = dict()
    num_group = 0
    with open(groupfile, 'r') as group:
        for line in group: # format: sp|A0A183|LCE6A_HUMAN   243797  2.40375
            group_info = line.strip("\n").split()
            seqID = group_info[0]
            if seqID[0] == 's':
                seqID = seqID.split('|')[1] # only the accession number (primary ID) is extracted, e.g. A0A183
                groupID = group_info[1] # groupID as integer, so that it can be ordered (by number, i.e. length)
                seq_dict[seqID] = {}
                # if groupID is not already a key in the dictionary, create
                if groupID not in group_dict:
                    group_dict[groupID] = [seqID];
                # otherwise, append the sequenceID to the list of sequences for this group
                else:
                    group_dict[groupID].append(seqID);
                num_group += 1
                if num_group % 10000 == 0:
                    sys.stderr.write("{} groups read so far\n".format(num_group))
            else:
                continue
    sys.stderr.write("\nnumber of groups is {}\n".format(len(group_dict)))
    return group_dict,seq_dict

def read_sprot_dat(sprot_dat_file, seq_dict):
    num_record = 0
    for record in SwissProt.parse(open(sprot_dat_file)): # Use Bio.SwissProt to parse the uniprot_sprot.dat file
        for seqID in record.accessions:
            if seqID in seq_dict:
                num_record += 1
                if num_record % 10000 == 0:
                    sys.stderr.write("{} records read so far\n".format(num_record))
                go_terms = [i[1][3:] for i in record.cross_references if i[0] == 'GO'] # GO terms ['GO:0031012', 'GO:0005576', 'GO:0004222', 'GO:0008270']
                organism = record.organism # organism name
                lineage = record.organism_classification # taxonomic classification ['Viruses', 'dsDNA viruses, no RNA stage', 'Iridoviridae', 'Chloriridovirus']
                tax_id = record.taxonomy_id[0] # taxonomy id '345201'
                gene_name, OLN, ORF = parse_GN(record.gene_name) # GN line,include gene names, ordered locus names, and ORF names 
                full_name, EC = parse_DE(record.description) # DE line with descriptive information. RecName, AltName (Full=, short=, EC=, ...)
                seq_dict[seqID] = {'organism' : organism, 'EC' : EC, 'gene_name' : gene_name,'OLN' : OLN, 'ORF' : ORF, 'GO': go_terms, 'KW': record.keywords, 'full_name': full_name, 'tax_id': tax_id, 'lineage': lineage} # map primary ID to annotation dictionary
            else:
                continue
    sys.stderr.write("\nnumber of sequences is {}\n".format(len(seq_dict)))
    return seq_dict

def parse_GN(string):
    ''' Parse GN (gene_name) record for a swissprot protein
        Return list (gene_name, OLN, ORF)
    '''
    gene_name = [x.lower() for x in get_value(string, "Name") + get_value(string, "Synonyms")]
    OLN = get_value(string, "OrderedLocusNames")
    ORF = get_value(string, "ORFNames")
    return gene_name, OLN, ORF

def parse_DE(string):
    ''' Parse DE (description) record for a swissprot protein
        Return list (full_name(including short_name), EC)
    '''
    full_name = get_value(string, "Full")
    short_name = get_value(string, "Short")
    EC = get_value(string, "EC")
    return full_name+short_name, EC

def get_value(string, field):
    ''' From a string in sprot.dat, find the value for a given field
        Return list.
        Format: field=val1, val2 {evidence}, ..., valn {evidence}; 
    '''
    values = []
    start_pos = 0
    while(True):
        field_pos = string.find(field+"=", start_pos)
        if field_pos == -1:
            break
        end_pos = string.find(';', field_pos) # ending position of the values
        values = values + strip_evidence(string[(field_pos+len(field) +1):end_pos]).split(',') # values are separated by comma
        start_pos = end_pos + 1
    return [x.strip() for x in values]

def strip_evidence(string):
    ''' Strip the {evidence} string from the value string
        Ex: val2 {evidence} -> val2
    '''
    return re.sub(r'{.*}', '',string.strip())
## =================================================================
## main function
## =================================================================
def main(argv=None):
    if argv is None:
        args = parser.parse_args()

    if args.verbose:
        sys.stderr.write( 'running verbosely\n')
        sys.stderr.write( 'input dat file is: {}\n'.format(args.dat_file))
        sys.stderr.write( 'input group file is: {}\n'.format(args.groupfile))
        sys.stderr.write( 'output file is: {}\n'.format(args.outputfile))

    ## display work start, and time record
    start_time = datetime.now()
    sys.stderr.write("\n===============================================================================\n")
    sys.stderr.write("Assigning annotations for each group\n")
    ## run program and generate the fasta files
    annotation(args.dat_file,args.groupfile,args.outputfile)
    ## display work end, and time record, elapsed time
    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write("\nTotal Elapsed Time = [%s] [seconds] \n" % duration)
    sys.stderr.write("===============================================================================\n")

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
