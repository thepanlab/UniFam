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

void read_group_sg(const string & file_name, unordered_map<string, UINT16> & seq_gp, map<UINT16, UINT16> & group_len_map){
	ifstream ingroup(file_name.c_str());
	if(!ingroup)
	{
		cerr << "Input group file " << file_name << " cannot be opened for reading." << endl;
		exit(1);
	}
	string seqID, field;
	UINT16 groupID, seq_len, group_len;
	seq_gp.clear();
	seq_gp.reserve(10000000);

	while(ingroup >> seqID)
	{
		ingroup >> groupID >> field;
		if(field.back() == '*'){        /* centroid sequence of the group */
			field.pop_back();
			group_len = stoul(field);
			seq_len = group_len;
			group_len_map[groupID] = group_len;
		}
		else
			seq_len = stoul(field);
		seq_gp[seqID] = seq_len;
	}
	cout << "Total number of groups: " << group_len_map.size() << endl;
	cout << "Total number of sequences: " << seq_gp.size() << endl;
	cout << "Hash load factor is " << seq_gp.load_factor() << endl;
}
