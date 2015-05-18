#ifndef READFILE_H
#define READFILE_H

/********************************************************************
 * ReadFile.h
 * Created by JJ Chai  on 05/13/2015.
 * Copyright (c) 2015 JJ Chai (ORNL). Allrights reserved.
********************************************************************/

#include <ctime>
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
#include <unordered_map>
#include <sys/stat.h>
#include "utils.hpp"

typedef unsigned int UINT16;
typedef unsigned long long int UINT64;
using namespace std;

/* Read group file (seqID groupID seqLen) and save information */
void read_group_sg(const string & file_name, unordered_map<string, UINT16> & seq_gp, map<UINT16, vector<UINT16> > & group_len_map);
void read_faidx(const string & file_name, unordered_map<string, UINT64> & seq_ofs);
void read_group(const string & group_file, unordered_map<string, vector<string> > & groupMap);

void filter_fasta(const string & fasta_file, const unordered_map<string, UINT16> & seq_gp, const string & out_fasta);
void get_seq(const string & fasta_file, const string & seqID, const UINT64 & offset, ofstream & fout);
void split_group(unordered_map<string, vector<string> > & groupMap, const size_t & size, const string & output_prefix);

void group_seqs(const string & group_name, const unordered_map<string, UINT64> & seq_ofs, const unordered_map<string, vector<string> > & groupMap,const string & fasta_file, const string & fasta_name, size_t & group_size);
int msa_hmm(const string & group_name, const size_t & group_size, const string & fasta_name, const string & hmm_file, const string & mafft_prefix, const string & hmm_prefix);
#endif //READFILE.HPP
