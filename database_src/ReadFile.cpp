/*
 * =====================================================================================
 *
 *       Filename:  ReadFile.cpp
 *
 *    Description:  Read different types of files and save them into memory with different 
 *    		    types of data structures
 *
 *        Version:  1.0
 *        Created:  05/13/2015 15:44:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  JJ Chai (jjchai), jjchai01@gmail.com
 *   Organization:  ORNL
 *
 * =====================================================================================
 */

#include "ReadFile.h"
using namespace std;
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  read_group_sg
 *  Description:  Read group file (seqID groupID seqLen) and save information: 
 *                hashmap: seqID -> groupID
 *                map: groupID -> length of centroid
 * =====================================================================================
 */

void read_group_sg(const string & file_name, unordered_map<string, UINT16> & seq_gp, map<UINT16, vector<UINT16> > & group_len_map){
	ifstream ingroup(file_name.c_str());
	if(!ingroup)
	{
		cerr << "Input group file " << file_name << " cannot be opened for reading." << endl;
		exit(1);
	}
	string seqID, field;
	UINT16 groupID, group_len, num_groups;
	seq_gp.clear();
	seq_gp.reserve(10000000);

	num_groups = 0;
	while(ingroup >> seqID)
	{
		ingroup >> groupID >> field;
		if(field.back() == '*'){        /* centroid sequence of the group */
			field.pop_back();
			group_len = stoul(field);
			group_len_map[groupID].push_back(group_len);
			group_len_map[groupID].push_back(num_groups);
			num_groups++;

		}
		seq_gp[seqID] = groupID;
	}
//	cout << "Total number of groups: " << group_len_map.size() << endl;
//	cout << "Total number of sequences: " << seq_gp.size() << endl;
//	cout << "Hash load factor is " << seq_gp.load_factor() << endl;
	ingroup.close();
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  filter_fasta
 *  Description:  Given a fasta file, and a sequence to group map, filter out the sequences
 *  	     	  that are not included in any group.
 * =====================================================================================
 */
void filter_fasta(const string & fasta_file, const unordered_map<string, UINT16> & seq_gp, const string & out_fasta)
{
	ifstream fin(fasta_file.c_str());
	if(!fin)
	{
		cerr << "Input group file " << fasta_file << " cannot be opened for reading." << endl;
		exit(1);
	}

	char first_char;
	string seqID, seq;
	bool print=false;
	ofstream fout(out_fasta.c_str()); /* output fasta file */
	while(!fin.eof()){
		fin.get(first_char);
		if(first_char == '>'){                  /* header line */
			print = false;
			getline(fin, seqID);
			fin.ignore(numeric_limits<streamsize>::max(),'\n');
			if(seq_gp.count(seqID) > 0){
				print = true;
				fout << ">" << seqID << endl;
			}
		}
		else if(print){
			getline(fin, seq);
			fout << seq << endl;
		}
		else{
			fin.ignore(numeric_limits<streamsize>::max(),'\n');
		}
	}
	fin.close();
}


void get_seq(const string & fasta_file, const string & seqID, const UINT64 & offset, ofstream & fout)
{
	ifstream fin(fasta_file.c_str());
	if(!fin)
	{
		cerr << "Input group file " << fasta_file << " cannot be opened for reading." << endl;
		exit(1);
	}
	fin.seekg(offset, fin.beg);
	fout << ">" << seqID << endl;
	char first_char;
	string seq;
	fin >> first_char;
	while(first_char != '>'){
		fout << first_char;
		getline(fin, seq);
		fout << seq << endl;
		fin >> first_char;
	}
	fin.close();
}


void read_faidx(const string & file_name, unordered_map<string, UINT64> & seq_ofs)
{
	ifstream fin(file_name.c_str());
	if(!fin)
	{
		cerr << "Input group file " << file_name << " cannot be opened for reading." << endl;
		exit(1);
	}
	seq_ofs.clear();
	seq_ofs.reserve(10000000);
	string seqID;
	UINT64 offset;
	UINT64 tmp;
	while(!fin.eof()){
		fin >> seqID >> tmp >> offset;
		fin.ignore(numeric_limits<streamsize>::max(),'\n');
		seq_ofs[seqID] = offset;
	}
//	cout << "Number of sequences: " << seq_ofs.size() << endl;
	fin.close();
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  read_group
 *  Description:  read group information into memory <groupID, seqIDs in the group>
 * =====================================================================================
 */
void read_group(const string & group_file, unordered_map<string, vector<string> > & groupMap)
{
	ifstream infile(group_file.c_str());
	if(!infile)
	{
		cerr << "Input group file " << group_file << " cannot be opened for reading." << endl;
		exit(1);
	}
	string seqID, group;
	groupMap.clear();

	while (infile >> seqID) {
		//infile >> seqID >> group;
		infile >> group;
		//cout << seqID << "\t" << group << endl;
		infile.ignore(numeric_limits<streamsize>::max(),'\n');
		groupMap[group].push_back(seqID);
		//cout << "push sequence " << seqID << endl;

	}
	infile.close();
}


void split_group(unordered_map<string, vector<string> > & groupMap, const size_t & size, const string & output_prefix)
{
	size_t num_group = groupMap.size();
	size_t g_index = 0;
	int piece = 0;
	unordered_map<string, vector<string> >::iterator it = groupMap.begin();
	while(it != groupMap.end()){
		string out_group = output_prefix + to_string(piece) + ".group";
		ofstream gout(out_group.c_str());
		size_t end = min(g_index + size, num_group);
		for(size_t i = g_index; i < end; i++){
			size_t group_size = (it->second).size();
			for(size_t j = 0; j < group_size; j++){
				gout << (it->second).at(j) << "\t" << it->first << endl;
			}
			it++;
			g_index++;
		}
		piece++;
		gout.close();
	}
}



void group_seqs(const string & group_name, const unordered_map<string, UINT64> & seq_ofs, const unordered_map<string, vector<string> > & groupMap,const string & fasta_file, const string & fasta_name, size_t & group_size)
{
	vector<string> seqIDs = groupMap.at(group_name); /* sequence IDs */
	group_size = seqIDs.size();
	ofstream fasta_out(fasta_name.c_str()); /* output fasta file stream */
	for(size_t i = 0; i < group_size; i++){
		string seqID = seqIDs.at(i);
		get_seq(fasta_file, seqID, seq_ofs.at(seqID), fasta_out);
	}
	fasta_out.close();
}


int msa_hmm(const string & group_name, const size_t & group_size, const string & fasta_name, const string & hmm_file, const string & mafft_prefix, const string & hmm_prefix)
{
	string mafft_cmd = "";
	string mafft_fasta = fasta_name;
	if(group_size > 1){
		mafft_fasta = fasta_name + ".mafft";
		mafft_cmd = mafft_prefix + " --auto --anysymbol --amino --quiet --thread 1 " + fasta_name + " > " + mafft_fasta;
//		cout << "mafft command: " << mafft_cmd << endl;
		const int res = system(mafft_cmd.c_str());
		if(res != 0)
		{
			cout << "mafft failed on " << group_name << endl;
			return res;
		}
		remove(fasta_name.c_str());
	}
	string hmm_cmd = hmm_prefix + " --amino --cpu 1 -o /dev/null -n " + group_name + " " + hmm_file + " " + mafft_fasta;
//	cout << "hmmbuild command: " << hmm_cmd << endl;

	const int res = system(hmm_cmd.c_str());
	if(res != 0)
		cout << "hmmbuild failed on " << group_name << endl;
	remove(mafft_fasta.c_str());
	return res;
}
