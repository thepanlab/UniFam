#ifndef PFCLUST_HPP
#define PFCLUST_HPP

/********************************************************************
 * pfclust.hpp
 *
 * pfclust.cpp includes functions used in the pipeline.
 * 
 * Created by JJ Chai  on 10/17/2013.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
********************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <map>
#include <sys/stat.h>
//#include <omp.h>
#include "directoryStructure.hpp"
#include "utils.hpp"
#include "runBioProg.hpp"
using namespace std;

namespace PfClust
{
#ifndef GROUP_SCORE
#define GROUP_SCORE
    // a structure to store a paired value <groupID, corresponding bit score>
    struct group_score {
        int groupID;
        float bitscore;
    };
#endif
//    // function to find the indices of the lower bound and upper bound
//    void find_index_vector(vector<int> & seqLens,int & ub, int & lb,int & index1,int & index2);
    
    // given a sequence length and the difference proportion, find the lower bound and upper bound
    void find_range(const int length, float prop);
    
    //given file name of hmm, do hmmsearch
    void do_hmmsearch(const string & hmmFilename, const string & outputDir,
                      map<string, string> & seqMap, const map<int, vector<string> > & lenMap,
                      const float & prop, const string & hmmsearch_cmd_prefix,
                      const bool & log, const bool & clean);
    
    // get all seqIDs from group_sg
    void get_seqIDs(const string & group_sg, map<string, int> & seqID_index,
                    const float & cutoff, int & num_seqs);
    
    //read fasta file and store sequences in a map <seqID => seq>
    // and the lengths of sequences in a vector of integers <length => seqIDs with this length>
    void read_fasta(const string & fastafile, map<string, string> & seqMap, map<int, vector<string> > & lenMap);
    
    //read fasta file and store sequences in a map <seqID => seq>
    void read_fastaseq(const string & fastafile, map<string, string> & seqMap);
    
    // read a (big) fasta file and write their length information to a file
    void get_fasta_len(const string & fastafile, const string &lengthfile);
    
    // read a (big) fasta file and store their length information in a vector
    void read_fasta_len(const string & fastafile, map<string, int> & seqLenMap);
    
    // parse the group file (with bit score attached)
    void parse_group(const string & filename, map<string, group_score> & scores, const bool & clean);
    
    //parse all the groups in the directory, and generate overall group file
    void parse_groupDir(const string & groupfileDir, const string & output, const bool & clean);
    
    //parse all the groups in the directory, and generate set covering input file
    void setGen_groupDir(const string & groupfileDir, const string & output, const string & groupIDfile,
                         const map<string, int> & seqID_index, const int & num_seqs, const float & cutoff,
                         const bool & clean);
    
    // parse the domain wise table from hmmsearch, find the indices and bit score of the target sequences
    // build map < seqID => bitwise score per residue >
    // no numeric indices, just seqID as id
    void parse_domtab(const string & filename, // domtab file name
                          map<string,float> & targetScores, const bool & clean); //resulting map
    // parse all domtab files in a directory
    void parse_domtabDir(const string & domtabDir, const string & groupfile, const bool & clean);
    void parse_domtabs(const vector<string> & domtabFilenames, const string & groupfile, const bool & clean);
    
    
    // parse domtab file generated from hmmsearch, with multiple hmms in the hmm database, instead of only one
    void parse_domtabfile(const string & domtabDir, const string & groupfile);
    // reorder group file by sequence ID
    void reorder_group(const string & filename, string & ordered_file);
    
    // generate group file by group from group file by sequence
    void sg_to_gs(const string & group_sg, const string & group_gs);
    
    // read group information into memory <groupID, seqIDs in the group>
    void read_group(const string & group_sg, map<int, vector<string> > & groupMap, vector<int> & groupNames);
    
    // read group file (with format seq_group), store the group info
    // void read_group(const string & group_sg, );
    
    // find the length of hmm from .hmm file
    int get_hmm_len(const string & hmmfile);
    
    // generate fasta files, one for each group
    void group_fasta(const map< int, vector<string> > & groups, const string & outputDir,
                     const string & singletonDir,
                     map<string, string> & seqMap);
    // generate fasta file for a single group
    void group_fasta_single(const int & groupID, const vector<string> & seqIDs,const string & nonSingletonDir,const string & singletonDir, map<string, string> & seqMap);
    
    // align sequences to create msa
    void align_fasta(const string & align_cmd_prefix, string & fastafile, string & msaDir, bool clean);
    
    // build hmm for every group. if singleton, build directly;
    // otherwise use mafft do alignment first, then build
    void hmmbuild_jj(const string & hmmbuild_cmd_firstpart,
                     const string & fastafile, bool singleton,
                     const string & hmmDir);
    
    // split big fasta file into small ones, depending on the number of nodes
    void split_fasta(const map<string, string> & seqMap, const map<int, vector<string> > & lenMap,
                     const map<int, vector<string> > & groupMap, const vector<int> & groupNames,
                     const int & num_nodes,  const string & outputDir);
    
    // parse blast6out file to get pairwise sequence identities
    void parse_blast6out(const string & filename, vector<float> & identities, const bool & clean);
}

#endif //PFCLUST.HPP
