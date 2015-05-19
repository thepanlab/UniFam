/*
 * =====================================================================================
 *
 *       Filename:  split_group.cpp
 *
 *    Description:  Split all the groups into chunks, for later processing (HMM building)
 *
 *        Version:  1.0
 *        Created:  05/18/2015 09:20:41 AM
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

void initializeArguments(int argc, char **argv, string & group_file, string & prefix, UINT64 & chunk_size);
void usage();
/* print usage */
void usage()
{
    cout << " [Usage]" << endl
    << "  split_group -g <group_file> -s <chunk size> -p <output prefix>" << endl
    << endl
    << " [Inputs]" << endl
    << "   -g <group_file> group file with format (seqID \\t groupID)" << endl
    << endl
    << " [Outputs]" << endl
    << "   -p <output prefix> prefix (including directory name) of the split group files" << endl
    << endl
    << " [Options]" << endl
    << "   -s <chunk size> number of groups in each chunk (default: 10000)" << endl
    << "   -f <fasta file> fasta with all sequences" << endl
    << "   -h/--help Display this help message" << endl
    << endl;
}

/*
 * Parse command line arguments
 * split_group -g <group_file> -s <size> -p <output prefix>
 */

void initializeArguments(int argc, char **argv,
                         string & group_file,
                         string & fasta_file,
			 string & prefix,
                         UINT64 & chunk_size)
{
    vector<string> Arguments;
    
    group_file  = "";
    fasta_file  = "";
    prefix = "";
    chunk_size = 10000;
    
    while(argc--)
        Arguments.push_back(*argv++);
    
    for(int i = 1; i <= (int)Arguments.size()-1; i++)
    {
        // input group file
        if (Arguments[i] == "-g") {
            group_file = Arguments[++i];
        }
        
	// input fasta file
	else if (Arguments[i] == "-f") {
            fasta_file = Arguments[++i];
        }
        
	// output group file prefix
	else if (Arguments[i] == "-p") {
            prefix = Arguments[++i];
        }
        
        // number of groups in each file
        else if (Arguments[i] == "-s")
            chunk_size = stoull(Arguments[++i]);
        
        // help, print Usage
        else if ((Arguments[i] == "-h") || (Arguments[i] == "--help"))
        {
            usage();
            exit(0);
        }
        else
        {
            cerr << "Unknown option " << Arguments[i] << endl << endl;
            usage();
            exit(1);
        }
    }
    // check required arguments
    if ((group_file == "") || (prefix == "")){
        cerr << "miss necessary parameter(s)"<<endl;
        cerr << "use -h/--help for help" << endl;
        exit(1);
    }
    cout << "File with all the sequences and their groups: " << group_file << endl;
    cout << "Fasta file with all the sequences: " << fasta_file << endl;
    cout << "Prefix for the output group file: " << prefix << endl;
    cout << "Number of groups in each output group file: " << chunk_size << endl;
    
}


int main(int argc, char ** argv){

	/********************** display work start and time record *********************/

	time_t start_time = time(0); /* starting time in epoch */
	string group_file, prefix, fasta_file="";
	UINT64 chunk_size;
	initializeArguments(argc, argv, group_file, fasta_file, prefix, chunk_size);
	cout << endl
		<< "============================================================================" << std::endl
		<< Utils::currentDateTime()
		<< endl;

	/********************** variable declaration *********************/

	unordered_map<string, vector<string> > groupMap;	

	/********************** do things *********************/
	cout << "Read group file " << group_file << endl;
	read_group(group_file, groupMap);
	if(fasta_file.size() > 0){
		cout << "Read fasta file " << fasta_file  << endl;
		unordered_map<string, string> seqMap;
		read_fastaseq(fasta_file, seqMap);
		cout << "Split group and fasta files.\n";
		split_group_fasta(groupMap, seqMap, chunk_size, prefix);
	}
	else{
		cout << "Split group file only.\n";
		split_group(groupMap, chunk_size, prefix);
	}

	/********************** display work end and time record *********************/

	time_t end_time = time(0); /* ending time in epoch */
	cout<< "Done!" << endl
		<< Utils::currentDateTime() << endl
		<< "Total elapsed time = " << double(end_time - start_time) << " [seconds]" << endl
		<< "============================================================================"
		<< endl << endl;

	return 0;
}
