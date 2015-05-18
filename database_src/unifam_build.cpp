/********************************************************
 * unifam_build.cpp
 
 * This is the final pipeline
 
 * Created by JJ Chai  on Thu May 14 16:03:00 EDT 2015. Last modified Thu May 14 EDT 2015
 * Copyright (c) 2015 JJ Chai (ORNL). Allrights reserved.
 
 ********************************************************/

#include "ReadFile.h"
#include "directoryStructure.hpp"
#include <algorithm>
#include <mpi.h>
#include <omp.h>

#define WORKTAG    1
#define DIETAG     2

using namespace std;
// slave process, everybody is slave...
void SlaveProcess(const string & fasta_file, const vector<string> & group_files, const string & tmpDir, const string &outputDir,
                  const string & mafft_cmd_path, const string & hmmbuild_cmd_path, 
		  const unordered_map<string, UINT64> & seq_ofs, const size_t & nThreads);
void MasterProcess(const size_t & num_groups);
/* print usage */
void usage()
{
    cout << " [Usage]" << endl
    << "  unifam_build [options] -f <fasta_file> -in <group_dir> -td <tmpDir> -od <outputDir> " << endl
    << "  -align <mafft_cmd_path> -build <hmmbuild_cmd_path>" << endl
    << endl
    << " [Inputs]" << endl
    << "   -f <fasta_file> fasta file including all the protein sequences" << endl
    << "   -in <group_dir> directory with all the group files" << endl
    << "   -align <mafft_cmd_path> full path for the alignment command (mafft)" << endl
    << "   -build <hmmbuild_cmd_path> full path for hmmbuild" << endl
    << endl
    << " [Outputs]" << endl
    << "   -td <tmpDir> directory for temporary intermediate files" << endl
    << "   -od <outputDir> directory for new group files" << endl
    << endl
    << " [Options]" << endl
    << "   -p <nThreads> number of threads for each process (default: all available threads)" << endl
    << "   -h/--help Display this help message" << endl
    << endl;
}

void initializeArguments(int argc, char **argv,
                         string & fasta_file, // fasta file with all sequences to extract
			 string & group_dir, // directory with all the group files
                         string & tmpDir, // temporary output directory
                         string & outputDir, // output directory for the group files
                         string & mafft_cmd_path,//
                         string & hmmbuild_cmd_path,
			 size_t & nThreads
)
{
    vector<string> Arguments;
    
    fasta_file  = "";
    group_dir  = "";
    tmpDir = "";
    outputDir = "";
    nThreads = 0;
    mafft_cmd_path = "mafft";
    hmmbuild_cmd_path = "hmmbuild";
    
    while(argc--)
        Arguments.push_back(*argv++);
    
    for(int i = 1; i <= (int)Arguments.size()-1; i++)
    {
        // input fastafile
        if (Arguments[i] == "-f") {
            fasta_file = Arguments[++i];
        }
        
	// group file directory
	else if (Arguments[i] == "-in") {
            group_dir = Arguments[++i];
        }
        
        // tmp outputDir
        else if (Arguments[i] == "-td")
            tmpDir = Arguments[++i];
        
        // group file outputDir
        else if (Arguments[i] == "-od")
            outputDir = Arguments[++i];
        
        // mafft
        else if (Arguments[i] == "-align")
            mafft_cmd_path = Arguments[++i];
        
        // hmmbuild
        else if (Arguments[i] == "-build")
            hmmbuild_cmd_path = Arguments[++i];
        
        // threads
        else if (Arguments[i] == "-p")
            nThreads = static_cast<size_t>(stoi(Arguments[++i]));
        
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
    if ((group_dir == "") || (fasta_file == "") || (tmpDir == "") || (outputDir == "") || (mafft_cmd_path == "")
        || (hmmbuild_cmd_path == ""))
    {
        cerr << "miss necessary parameter(s)"<<endl;
        cerr << "use -h/--help for help" << endl;
        exit(1);
    }
    if (nThreads == 0)
	    nThreads = omp_get_max_threads();

//    cout << "Input fasta file is: " << fasta_file << endl;
//    cout << "Input directory with group files is: " << group_dir << endl;
//    cout << "Temporary directory is: " << tmpDir << endl;
//    cout << "Output directory is: " << outputDir << endl;
//    cout << "mafft is at: " << mafft_cmd_path << endl;
//    cout << "hmmbuild is at: " << hmmbuild_cmd_path << endl;
//    cout << "Number of threads for each proecess is: " << nThreads << endl;
    
}


// main function

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  main function
 * =====================================================================================
 */
int main(int argc, char ** argv){
    
    /********************** variable declaration *********************/
    
    string fasta_file, group_dir, mafft_cmd_path, hmmbuild_cmd_path, tmpDir, outputDir;
    vector<string> group_files;
    unordered_map<string, UINT64> seq_ofs; 
    size_t nThreads; 
    int myid, num_nodes;
    double start_time, finish_time;
    
    /********************** initialize MPI  *********************/
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
    num_nodes = num_nodes - 1;
    /********************** display work start and time record *********************/
    if (myid == 0) {
        start_time = MPI_Wtime();
        
        cout << endl
        << "============================================================================"
        << endl << Utils::currentDateTime() << endl
        << " Begin pipeline" << endl
        << " [Step 1] Initialize commandline arguments -> " << endl;
    }
    
    /************ initialize arguments, for input, output, and options ***********/
    
//    cout << "Maximum threads possible: " << nThreads << endl;
    initializeArguments(argc, argv, fasta_file, group_dir, tmpDir, outputDir, mafft_cmd_path, hmmbuild_cmd_path, nThreads);
    Utils::getFiles(group_dir, "group", group_files);
    size_t num_groups = group_files.size();
    
    Utils::mkdirIfNonExist(tmpDir);
    Utils::mkdirIfNonExist(outputDir);

    read_faidx(fasta_file+".fai", seq_ofs);
    
    /* echo the arguments */
    if (myid == 0) {
        cout<< "   Fasta file with all proteins is : " << fasta_file << ", total sequences: " << seq_ofs.size() << endl
            << "   Directory with all the group files is : " << group_dir << ", which has " << num_groups << " groups" << endl
            << "   Temporary output directory is : " << tmpDir << endl
            << "   Output directory is : " << outputDir << endl
            << "   Number of threads per process: " << nThreads << endl
            << "   Full path for the multiple sequence alignment command is : " << mafft_cmd_path << endl
            << "   Full path for hmmbuild command is : " << hmmbuild_cmd_path << endl << flush;
        cout<< " [Step 1] Done!" << endl << endl;
    }
    
    /********************** generate data for each node, and find chunk information *********************/
    
    if(num_nodes > 0){
	    if (myid == 0) {
		    MasterProcess(num_groups);
		    cout << "  Running with " << num_nodes << " slave processes" << endl;
	    }

	    else
		    SlaveProcess(fasta_file, group_files, tmpDir, outputDir, mafft_cmd_path, hmmbuild_cmd_path, seq_ofs, nThreads);
    }
    else{ /* If there is only one node, running in multi-thread mode */
	    cout << "Only 1 process, running with multiple threads. \n";
	    for(size_t myGroup = 0; myGroup < num_groups; myGroup++){
		    string group_file = group_files.at(myGroup);
		    unordered_map<string, vector<string> > groupMap;	
		    read_group(group_file, groupMap);
		    vector<string> group_names;
		    for(unordered_map<string, vector<string> >::iterator it = groupMap.begin(); it != groupMap.end(); it++)
			    group_names.push_back(it->first);

		    omp_set_num_threads(nThreads);
		    string hmm_dir = outputDir + Utils::getPathSeparator() + to_string(myGroup);
		    Utils::mkdirIfNonExist(hmm_dir);
		#pragma omp parallel for
		    for (size_t i = 0; i < group_names.size(); i++) {
			    string group_name = group_names.at(i);
			    int res = 0;
			    size_t group_size;
			    string tmp_fasta = tmpDir + Utils::getPathSeparator() + group_name + ".fasta";
			    string hmm_name = hmm_dir + Utils::getPathSeparator() + group_name + ".hmm";
			    group_seqs(group_name, seq_ofs, groupMap, fasta_file, tmp_fasta, group_size);
			    res = msa_hmm(group_name, group_size, tmp_fasta, hmm_name, mafft_cmd_path, hmmbuild_cmd_path);

		    }
	    }
    }
    
    /********************** display work end and time record *********************/
    
    if( myid == 0 )
    {
        cout << " Done!" << endl;
        finish_time = MPI_Wtime();
        
        // display work end and time record
        cout << Utils::currentDateTime() << " all nodes done "<< endl
        << "Total Elapsed Time =  "
        << double(finish_time - start_time) << " [seconds]" << endl
        << "============================================================================"
        << std::endl << std::endl;
    }
    
     /*************************** Finalize MPI  *******************************/
    MPI_Finalize();
    return 0;
}

void MasterProcess(const size_t & num_groups)
{
    int result, total_proc;
    size_t num_slave, num_MPI, currentGroup;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &total_proc);
    
    num_slave = total_proc - 1;
    num_MPI = (num_groups <= num_slave) ? num_groups : num_slave; /* number of slave MPI processes to start */
    
    for (size_t i = 1; i <= num_MPI; i++) {
        currentGroup = i - 1;
        MPI_Send(&currentGroup, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
        //cout << "   " << "0 : send signal to " << i << endl;
    }
    
    if (num_groups > num_slave) {
        currentGroup = num_slave;
        
        while (currentGroup < num_groups) {
            MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            //cout << "   " << "0: receive signal from " << status.MPI_SOURCE << endl;
            MPI_Send(&currentGroup, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            currentGroup++;
        }
    }
    // all zeros for dying command
    int die = 0;
    for (size_t i = 1 ; i <= num_slave; i++) {
        MPI_Send(&die, 1, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
    }
}



void SlaveProcess(const string & fasta_file, const vector<string> & group_files, const string & tmpDir, const string &outputDir,
                  const string & mafft_cmd_path, const string & hmmbuild_cmd_path, 
		  const unordered_map<string, UINT64> & seq_ofs, const size_t & nThreads)
{
    MPI_Status status;
    int myid;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    while (true) {
        int myGroup;
        MPI_Recv(&myGroup, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == DIETAG) {
            cout << "   " << myid << ": I'm going to die " << endl;
            break;
        }
        
        /* do the pipeline process for all the groups in this chunk */
	string group_file = group_files.at(myGroup);
	unordered_map<string, vector<string> > groupMap;	
	read_group(group_file, groupMap);
	vector<string> group_names;
	for(unordered_map<string, vector<string> >::iterator it = groupMap.begin(); it != groupMap.end(); it++)
		group_names.push_back(it->first);
        
	string hmm_dir = outputDir + Utils::getPathSeparator() + to_string(myGroup);
	Utils::mkdirIfNonExist(hmm_dir);
	omp_set_num_threads(nThreads);
        #pragma omp parallel for
        for (size_t i = 0; i < group_names.size(); i++) {
            string group_name = group_names.at(i);
            int res = 0;
	    size_t group_size;
            string tmp_fasta = tmpDir + Utils::getPathSeparator() + group_name + ".fasta";
            string hmm_name = hmm_dir + Utils::getPathSeparator() + group_name + ".hmm";
	    group_seqs(group_name, seq_ofs, groupMap, fasta_file, tmp_fasta, group_size);
	    res = msa_hmm(group_name, group_size, tmp_fasta, hmm_name, mafft_cmd_path, hmmbuild_cmd_path);
            
        }
        /************************ parse domtab files *********************************/
        cout << "   " << myid << ": all hmms generated for group file " << group_file << endl;
        
        MPI_Send(0,0,MPI_INT,0,0,MPI_COMM_WORLD);
    }
}
