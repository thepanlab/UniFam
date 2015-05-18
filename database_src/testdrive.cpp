/*
 * =====================================================================================
 *
 *       Filename:  testdrive.cpp
 *
 *    Description:  Test new functions and classes
 *
 *        Version:  1.0
 *        Created:  05/13/2015 16:03:27
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
int main(int argc, char ** argv){

	/********************** display work start and time record *********************/

	time_t start_time = time(0); /* starting time in epoch */
	cout << endl
		<< "============================================================================" << std::endl
		<< Utils::currentDateTime()
		<< endl
		<< "Beginning  -> " << endl;

	/********************** variable declaration *********************/

	string fasta_file = argv[1];
	string group_file = argv[2];
	string group_name = argv[3];
//	size_t size = stoi(argv[3]);
//	string seqID = argv[2];
//	UINT16 offset = stoull(argv[3]);
	unordered_map<string, vector<string> > groupMap;	
	unordered_map<string, UINT64> seq_gp; 
	size_t group_size;
//	map<UINT16, vector<UINT16> > group_len_map;

	/********************** test run of function *********************/

//	filter_fasta(fastafile, seq_gp, outputfasta);
//	get_seq(fasta_file, seqID, offset);
	cout << "Read faidx file " << fasta_file + ".fai" << endl;
	read_faidx(fasta_file+".fai", seq_gp);
	cout << "Read group file " << group_file << endl;
	read_group(group_file, groupMap);
//	split_group(groupMap, size, prefix);
	cout << "Generate fasta for group " << group_name << endl;
	group_seqs(group_name, seq_gp, groupMap, fasta_file, "tmp_fasta", group_size);
	cout << "Group size is " << group_size << endl;
	cout << "mafft and hmmbuild" << endl;
	int res = msa_hmm(group_name, group_size, "tmp_fasta", "tmp.hmm", "mafft", "hmmbuild");
	/********************** display work end and time record *********************/

	time_t end_time = time(0); /* ending time in epoch */
	cout<< "Done!" << endl
		<< Utils::currentDateTime() << endl
		<< "Total elapsed time = " << double(end_time - start_time) << " [seconds]" << endl
		<< "============================================================================"
		<< endl << endl;


	return 0;
}
