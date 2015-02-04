#!/usr/bin/python

'''
gbkGen.py

Generate Genbank file from other information

Created by JJ Chai on Wed Jan 21 14:19:17 EST 2015
Last modified 01/21/2015
Copyright (c) 2015 JJ Chai (ORNL). All rights reserved.

'''
# Import Python modules
import ConfigParser
import argparse
import sys, os, re
from datetime import datetime
import time

## Version
version_str = "0.1"
''' 0.1.  first version
'''


## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="Generate GenBank file from necessary input files (written for Sagar)",
                                 prog = 'gbkGen', #program name
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

## input files 
parser.add_argument("-gtf",help="genemark gtf output",dest='gtfFile',required=True)
parser.add_argument("-faa",help="fasta file with annotated proteins",dest='faaFile',required=True)
parser.add_argument("-fna",help="fasta file with annotated proteins",dest='fnaFile',required=True)
parser.add_argument("-annot",help="annotation flat file for the proteins",dest='annotFile',required=False) # If proteins in faaFile are not annotated, user can supply this flat annotation file instead


## output file, determined from the prefix of the input file by default if not specified
parser.add_argument("-o",help="output gbk file's prefix",dest='outputPrefix',required=False)

## =================================================================
## gbkGen function
## =================================================================
def gbkGen(gtfFile, faaFile, fnaFile, outputPrefix, annotFile=None):
    # read gtf file output from genemark, and save the gene information in a dictionary
    previous_contig = "null"
    contigs_seq = read_fna(fnaFile)
    # read protein sequences and their annotation
    if annotFile is not None: # annotation file is provided as flat file in addition to the faa file
        annot = read_annot(annotFile)
    else: # annotation is attached to the header line of the protein sequences
        annot = dict()
        sys.stderr.write("Warning: no annotation file provided, hence no genes and gene product available.\n")
    proteins_seq = read_fna(faaFile)
    datenow = time.strftime("%d-%b-%Y")
    #print datenow
    with open(gtfFile) as gtf:
        for gene_record in read_gtf(gtf):
            contig,location_string, geneID = extract_info(gene_record)
            if contig != previous_contig:
                # close the file for the old contig
                if previous_contig != "null":
                    gbk.write('{0:6}\n'.format("ORIGIN")) # write the whole genome sequence in the end
                    #gbk.write('{0:>9} {}\n'.format(1, 'acgt...'))
                    write_fna(contigs_seq[contig], gbk) # write genome nucleotide sequence at the end
                    gbk.write("//\n")
                    gbk.close() 
                    #break

                previous_contig = contig # current contig becomes "previous_contig"

                # write to a new file for a new contig
                gbk = open(outputPrefix + "_" + contig + ".gbk", "w") 
                gbk.write("{:12}{}\t{}\t{}\t{}\t{}\n".format("LOCUS", "SHORTNAME", "XX bp", "DNA", "PLN", datenow))
                gbk.write("{:12}{}\n".format("DEFINITION","definition, details about the contig/genome"))
                gbk.write("{:12}{}\n".format("ACCESSION", "unknown"))
                gbk.write("{:12}{}\n".format("KEYWORDS", "."))
                gbk.write("FEATURES{0:13}Location/Qualifiers\n".format(' ')) # FEATURES             Location/Qualifiers

                gbk.write("{:5}{:16}1..{}\n".format(' ','source', "length")) # several placeholders
                gbk.write("{:21}{}\n".format(' ', '/mol_type="genomic DNA"')) 
                gbk.write("{:21}{}\n".format(' ', '/organism="organism name"')) 
                gbk.write("{:21}{}\n".format(' ', '/note="note"')) 

                locus_tag = 1 # reset the locus tag to 1
            gbk.write("{:5}{:16}{}\n".format(' ','CDS', location_string)) # CDS line
            gbk.write('{:21}{}"{}{:0>5}"\n'.format(' ','/locus_tag=','LOCUS_', locus_tag )) 
            protein_seq = '/translation="' + proteins_seq[geneID] + '"'
            write_protein_seq(protein_seq, gbk)
            gbk.write('{:21}{}{}\n'.format(' ','/transl_table=','1')) 
            write_gbk(geneID, annot, gbk)
            locus_tag = locus_tag + 1
        gbk.write('{0:6}\n'.format("ORIGIN")) # write the whole genome sequence in the end
        #gbk.write('{0:>9} {}\n'.format(1, 'acgt...'))
        write_fna(contigs_seq[contig], gbk) # write genome nucleotide sequence at the end
        gbk.write("//\n")
        gbk.close() 



def write_protein_seq(seq, gbk):
    ''' write translation of nucleotide (into amino acids) to the gbk file, with 21 spaces up front and 58 chars every line '''
    cols = 58
    i = 0
    len_seq = len(seq)
    while i < len(seq):
        gbk.write('{:21}{}\n'.format(' ', seq[i : min(len_seq, i + cols)]))
        i = i + cols

## =================================================================
## write gbk file for a contig/chromosome/plasmid,
## protein coding gene part
## =================================================================
def write_gbk(geneID, annot, gbk):
    ''' Write a contig's protein annotation in the corresponding file
        Input:
            1. annot_prot: annotation for all the proteins in this contig
            2. gbk: file object to write
    '''
    ## ## /gene="gap",(required) /alt_name="gad", /alt_name="foo"
    if geneID in annot:
        annot_prot = annot[geneID]
        gene_names = annot_prot["gene_name"].split(":")
        if gene_names[0] != "NA":
            gbk.write('{0:21}{1:}\n'.format(' ','/gene="' + gene_names[0] + '"'))
            for k in xrange(1,len(gene_names)):
                gbk.write('{0:21}{1:}\n'.format(' ','/alt_name="' + gene_names[k] + '"'))
        else: # if no gene name is available, write this line
            gbk.write('{0:21}{1:}\n'.format(' ','/gene=""'))
        
        ## /EC_number="" (recommended)
        ECnames = annot_prot["ECname"].split(":")
        if ECnames[0] != "NA":
            for k in xrange(0,len(ECnames)):
                gbk.write('{0:21}{1:}\n'.format(' ','/EC_number="' + ECnames[k].strip("EC ") + '"'))
        
        ## /product="" (required) /product_comment="" (optional)
        product_names = annot_prot["full_name"].split(":")
        gbk.write('{0:21}{1:}\n'.format(' ','/product="' + product_names[0] + '"'))
        for k in xrange(1,len(product_names)):
            gbk.write('{0:21}{1:}\n'.format(' ','/product_comment="' + product_names[k] + '"'))
        
        ## /db_xref="GO:xxxxxx" (optional)
        GOs = annot_prot["GO"].split(":")
        if GOs[0] != "NA":
            for k in xrange(len(GOs)):
                gbk.write('{0:21}{1:}\n'.format(' ','/db_xref="GO:' + GOs[k] + '"'))

        ## /gene_comment="" (optional)
        ## include information from fields OLN, ORF, KW
        gene_comment = annot_prot["OLN"] + ":" + annot_prot["ORF"] + ":" + annot_prot["KW"]
        gene_comment = re.sub("NA:","",gene_comment)
        gene_comment = re.sub(":?NA","",gene_comment) # remove trailing NA, or :NA
        if gene_comment != "":
            gene_comment = re.sub(":",";\n" + " "*21, gene_comment)
            gbk.write('{0:21}{1:}\n'.format(' ','/note="' + gene_comment + '"'))
    else:
        gbk.write('{0:21}{1:}"{2:}"\n'.format(' ','/gene=',geneID))
        gbk.write('{0:21}{1:}\n'.format(' ','/product="' + 'hypothetical protein"'))
        gbk.write('{:21}{}"{}"\n'.format(' ','/note=','not annotated by UniFam')) # CDS line

    

def read_annoted_faa(faaFile):
    ''' Read fasta format file including proteins, and has annotations in the header line.
        Input:  faaFile - as described in the line above
        Output: proteins_seqs - dictionary: proteinID -> sequence
                annot - dictionary: proteinID -> annotation information ( A little redundant but not too bad if the files are not that big )
    '''
    return

## =================================================================
## function to read annotation file for a database
## ================================================================
def read_annot(annot_file):
    ''' Read the annotation file, save in a dictionary and return.
        
        Input: annotation file, with fields: groupID, org, ECname, gene_name, OLN, ORF, full_name, GO, KW
        
        Output: dictionary object
    '''
    annot_dict = dict()
    with open(annot_file,'r') as annot:
        for line in annot:
            line = line.strip("\n").split("\t") # strip EOF
            if line[0] == "geneID":  # skip header line
                next(annot) 
            if len(line) > 7:
                annot_dict[line[0]] = {'org' : line[2], 'ECname' : line[3], 'gene_name' : line[4], 'OLN' : line[5], 'ORF' : line[6],\
                'full_name' : line[7], 'GO' : line[8], 'KW' : line[9]}
            else:
                annot_dict[line[0]] = {'ECname' : line[2], 'gene_name' : line[1], 'full_name' : line[3], 'GO' : line[5], 'KW' : line[4]}
    return annot_dict


## =================================================================
## function to read the fna file of contig sequences into dictionary
## ================================================================
def read_fna(fna_file,trim=True,reverse=False):
    ''' Read .fna file with all DNA sequences, and return dictionary object: contigID => contigSeq
        Input:
        1. fna_file: input fna file
        2. trim: trim ID at first space
        3. reverse: if true, seq => ID, instead of ID => seq
        Output
        1. dictionary object
    '''
    fna = open(fna_file,'r')
    ## all the contigs in a list
    contigs = fna.read().split('\n>')
    fna.close()
    
    contigs[0] = contigs[0][1:] # trim '>' for the first contig
    # input contigs info: contigID => sequence for this contig
    input_info = dict()
    for contig in contigs:
        if trim:
            contigID = contig.split()[0] # ID of this contig
        else:
            contigID = contig.split("\n",1)[0] # full header
        contigSeq = contig.split("\n",1)[1] # its corresponding sequence
        contigSeq = re.sub("\r","",contigSeq) # remove windows style EOF
        contigSeq = re.sub("\n","",contigSeq) # remove unix style EOF

        if reverse: # contigSeq => contigID
            # what if two proteins have the same sequences??
            # store the corresponding IDs in a list
            if input_info.has_key(contigSeq):
                input_info[contigSeq].append(contigID)
            else:
                input_info[contigSeq] = [contigID]
        else:
            input_info[contigID] = contigSeq # link them together in the dictionary
    return input_info

## =================================================================
## write the nucleotide sequence to file
## =================================================================
def write_fna(seq, gbk):
    i = 0
    len_seq = len(seq)
    while i < len_seq:
        seq_60bp = seq[i:min(len_seq, i+60)].lower()
        j = 0
        gbk.write('{:>9}'.format(i+1))
        len_subseq = len(seq_60bp)
        while j < len_subseq: # print 60 bps on a line
            gbk.write(' {}'.format(seq_60bp[j:min(len_subseq, j+10)]))
            j = j + 10
        gbk.write("\n")
        i = i + 60

## =================================================================
## read_gtf function
## Seems like .gtf files have the contigs sorted, not .gff files
## =================================================================
def read_gtf(gtf):
    ''' Read the .gtf file generated by genmark_hmm and find the gene names and the coding region.
        Use yield function to generate a gene record a time, memory print is samller this way.
        Input:  gtf - file object opening .gtf file
        Output: a list [] with information for a gene
    '''
    CDS = False
    record = ""
    line = gtf.readline() # read a line from the gtf file
    while line != '':
        print line

        # ignore lines starting with # or empty lines
        if line[0] == "#" or line=="\n": # ignore comment lines or blank lines
            line = gtf.readline()
            continue

        feature = line.strip("\n").split("\t")[2] # each field of the line
        record = record + line
        if not CDS and feature == 'CDS':
            CDS = True
        elif CDS and feature != 'CDS':
            yield record[:-1] # do no include the last \n
            record = ""
            CDS = False
        line = gtf.readline()
            
def extract_info(gtf_record):
    start = False
    stop = False
    CDS_start = 1
    locations = []
    lines = gtf_record.split('\n') # split the record by EOL
    line = lines[0]
    contig, source, feature, start_pos, end_pos, score, strand, frame, attribute = line.split("\t") # each field of the line
    geneID = attribute[(attribute.find('"') + 1):attribute.find('"',attribute.find('"')+1)] # gene id
    if feature == "start_codon":
        start = True
        start_codon_pos = start_pos
    elif feature == "CDS":
        CDS_start = 0

    line = lines[-1]
    contig, source, feature, start_pos, end_pos, score, strand, frame, attribute = line.split("\t") # each field of the line
    if feature == "stop_codon":
        stop = True
        stop_codon_pos = end_pos

    for i in xrange(CDS_start,len(lines)-1): # CDS line
        line = lines[i]
        contig, source, feature, start_pos, end_pos, score, strand, frame, attribute = line.split("\t") # each field of the line
        if feature != "CDS": 
            sys.stderr.out("formatting problem:\n{}\n".format(gtf_record)) # all these lines should have feature "CDS"
            break
        location_string = start_pos + ".."
        if not start or start_pos != start_codon_pos or frame != '0': # partial at the start
            location_string = "<" + location_string
        if i == len(lines)-2: # If this line is the last line of the CDS lines
            if not stop or end_pos != stop_codon_pos or frame != '0':
                location_string = location_string + ">" + end_pos
            else:
                location_string = location_string + end_pos
        else: # if not, then treat it as partial on the end
            location_string = location_string + ">" + end_pos
        # check the strand and see if it's on the complement strand
        if strand == "-":
            location_string = "complement(" + location_string + ")"
        locations.append(location_string)
    print locations
    return [contig, "join(" + ','.join(locations) + ")" if len(locations) > 1 else locations[0], geneID]


## =================================================================
## main function
## =================================================================
def main(argv=None):
    if argv is None:
        args = parser.parse_args()
    
    if args.outputPrefix == None:
        args.outputPrefix = os.path.abspath(args.faaFile).split('.')[0]
    ## print some information
    if args.verbose:
        sys.stdout.write('running verbosely\n')
        sys.stdout.write('genemark gtf file is: {}\n'.format(args.gtfFile))
        sys.stdout.write('protein fasta file is: {}\n'.format(args.faaFile))
        if args.annotFile != "":
            sys.stdout.write('annotation file is: {}\n'.format(args.annotFile))
        sys.stdout.write('nucleotide fasta file is: {}\n'.format(args.fnaFile))
        sys.stdout.write('output file prefix is: {}\n'.format(args.outputPrefix))
    else:
        sys.stdout.write('\n')

    # display work start, and time record
    start_time = datetime.now()
    sys.stderr.write("\n===============================================================================\n")
    sys.stderr.write("Start running: \n")

    # Annotating with UniFam
    gbkGen(args.gtfFile, args.faaFile, args.fnaFile, args.outputPrefix, args.annotFile)

    ## display work end, and time record, elapsed time
    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write("\nTotal Elapsed Time = [%s] \n" % duration)
    sys.stderr.write("===============================================================================\n")

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
