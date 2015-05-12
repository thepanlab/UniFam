#include "utils.hpp"
#include "pfclust.hpp"
#include "runBioProg.hpp"
#include <algorithm>
#include <mpi.h>
#include <omp.h>
/********************************************************
 * pipeline.cpp
 
 * This is the final pipeline
 
 * Created by JJ Chai  on 10/27/2013. Last modified 11/13/2013.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
 
 ********************************************************/

#define WORKTAG    1
#define DIETAG     2

using namespace std;
// slave process, everybody is slave...
void SlaveProcess(const string & inputDir, const string & tmpDir, const string & outputDir,
                  const string & mafft_cmd_path, const string & hmmbuild_cmd_path);
void MasterProcess(const vector<vector<int> > & chunk_info);
/* print usage */
void usage()
{
    cout << " [Usage]" << endl
    << "  pipeline [options] -in <inputDir> -td <tmpDir> -gd <outputDir> " << endl
    << "  -align <mafft_cmd_path> -build <hmmbuild_cmd_path>" << endl
    << endl
    << " [Inputs]" << endl
    << "   -in <inputDir> directory that includes all input files" << endl
    << "   -align <mafft_cmd_path> full path for the alignment command (mafft)" << endl
    << "   -build <hmmbuild_cmd_path> full path for hmmbuild" << endl
    << endl
    << " [Outputs]" << endl
    << "   -td <tmpDir> directory for temporary intermediate files" << endl
    << "   -gd <outputDir> directory for new group files" << endl
    << endl
    << " [Options]" << endl
    << "   -h/--help Display this help message" << endl
    << endl;
}

/*
 * Parse command line arguments
 * pipeline [options] -f <fastafile> -g <groupfile> -td <tmpDir> -align <mafft_cmd_path>
 */

void initializeArguments(int argc, char **argv,
                         string & inputDir, // fasta file with all sequences to extract
                         string & tmpDir, // temporary output directory
                         string & outputDir, // output directory for the group files
                         string & mafft_cmd_path,//
                         string & hmmbuild_cmd_path
)
{
    vector<string> Arguments;
    
    inputDir  = "";
    tmpDir = "";
    outputDir = "";
    mafft_cmd_path = "";
    hmmbuild_cmd_path = "";
    
    while(argc--)
        Arguments.push_back(*argv++);
    
    for(int i = 1; i <= (int)Arguments.size()-1; i++)
    {
        // input fastafile
        if (Arguments[i] == "-in") {
            inputDir = Arguments[++i];
        }
        
        // tmp outputDir
        else if (Arguments[i] == "-td")
            tmpDir = Arguments[++i];
        
        // group file outputDir
        else if (Arguments[i] == "-gd")
            outputDir = Arguments[++i];
        
        // mafft
        else if (Arguments[i] == "-align")
            mafft_cmd_path = Arguments[++i];
        
        // hmmbuild
        else if (Arguments[i] == "-build")
            hmmbuild_cmd_path = Arguments[++i];
        
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
    if ((inputDir == "") || (tmpDir == "") || (outputDir == "") || (mafft_cmd_path == "")
        || (hmmbuild_cmd_path == ""))
    {
        cerr << "miss necessary parameter(s)"<<endl;
        cerr << "use -h/--help for help" << endl;
        exit(1);
    }
    
}


// main function

int main(int argc, char ** argv){
    
    /********************** variable declaration *********************/
    
    string fastafile; //input files
    string groupfile;
    string mafft_cmd_path, hmmbuild_cmd_path;
    string inputDir, tmpDir, outputDir;

    
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
    
    initializeArguments(argc, argv, inputDir, tmpDir,outputDir, mafft_cmd_path, hmmbuild_cmd_path);
    
    Utils::mkdirIfNonExist(tmpDir);
    
    Utils::mkdirIfNonExist(outputDir);
    
    /* echo the arguments */
    if (myid == 0) {
        cout<< "   Input directory is : " << inputDir << endl
            << "   Temporary output directory is : " << tmpDir << endl
            << "   Output directory is : " << outputDir << endl
            << "   Full path for the MSA command is : " << mafft_cmd_path << endl
            << "   Full path for hmmbuild command is : " << hmmbuild_cmd_path << endl << flush;
        cout<< " [Step 1] Done!" << endl << endl
            << " [Step 2] Get the information on small chunks -> " << endl;
    }
    
    /********************** generate data for each node, and find chunk information *********************/
    
    if (myid == 0) {
        /* NOTE: cannot make vector of arrays!!! */
        vector<vector<int> > chunk_info;
        // file that has information on all the chunks
        ifstream mychunk((inputDir+Utils::getPathSeparator()+"chunk.info").c_str());
        
        int big_chunk, chunk_start, chunk_end;
        
        while (mychunk >> big_chunk) {
            mychunk >> chunk_start >> chunk_end;
            vector<int> this_chunk;
            this_chunk.push_back(big_chunk);
            this_chunk.push_back(chunk_start);
            this_chunk.push_back(chunk_end);
            chunk_info.push_back(this_chunk);
        }
        mychunk.close();
        cout << "   number of chunks is: " << chunk_info.size() << endl;
        cout<< " [Step 2] Done!" << endl << endl
        << " [Step 3] Slaves work on each small chunk -> " << endl;
        
        // echo the input, double check if the file was read correctly
        for (vector<vector<int> >::iterator it = chunk_info.begin(); it != chunk_info.end(); it++)
            cout << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << endl;
        
        MasterProcess(chunk_info);
    }
    
    else
        SlaveProcess(inputDir, tmpDir, outputDir, mafft_cmd_path, hmmbuild_cmd_path);
    
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

void MasterProcess(const vector<vector<int> > & chunk_info)
{
    int i, num_MPI, currentChunk, total_proc, num_slave, result;
    int chunk_num = (int) chunk_info.size(); /* total number of small chunks */
    cout << "   number of chunks is: " << chunk_num << endl;
    
    int chunk_send[4]; /* integer array to send to slave processes, cannot send vector with unknown length !!! */
    
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &total_proc);
    
    num_slave = total_proc - 1;
    num_MPI = (chunk_num <= num_slave) ? chunk_num : num_slave; /* number of slave MPI processes to start */
    
    for (i = 1; i <= num_MPI; i++) {
        currentChunk = i - 1;
        
        vector<int> this_chunk_info = chunk_info.at(currentChunk); /* big_chunk_num start end */
        this_chunk_info.push_back(currentChunk); /* index of this small chunk */
        
        /* assemble the int array to be sent */
        for (int j = 0 ; j < 4; j++) {
            chunk_send[j] = this_chunk_info.at(j);
        }
        
        /* send the array to a slave proc: <chunk, chunk_gp_start, chunk_gp_end, sub_chunk_ind> */
        MPI_Send(&chunk_send, 4, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
        //cout << "   " << "0 : send signal to " << i << endl;
    }
    
    if (chunk_num > num_slave) {
        currentChunk = num_slave;
        
        while (currentChunk < chunk_num) {
            MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            vector<int> this_chunk_info;
            this_chunk_info = chunk_info.at(currentChunk);
            this_chunk_info.push_back(currentChunk);
            for (int j = 0 ; j < 4; j++) {
                chunk_send[j] = this_chunk_info.at(j);
            }
            //cout << "   " << "0: receive signal from " << status.MPI_SOURCE << endl;
            MPI_Send(&chunk_send, 4, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            currentChunk++;
        }
    }
    // all zeros for dying command
    for (i = 0 ; i < 4; i++) {
        chunk_send[i] = 0;
    }
    for (i = 1 ; i <= num_slave; i++) {
        MPI_Send(&chunk_send, 4, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
    }
}



void SlaveProcess(const string & inputDir, const string & tmpDir, const string &outputDir,
                  const string & mafft_cmd_path, const string & hmmbuild_cmd_path)
{
    MPI_Status status;
    int myid;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    bool clean = true; /* clean up input file after output is generated? */
    
    while (true) {
        int this_chunk_info[4];
        MPI_Recv(&this_chunk_info, 4, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == DIETAG) {
            cout << "   " << myid << ": I'm going to die " << endl;
            //string newgroup = outputDir + Utils::getPathSeparator() + "newgroup." + Utils::intToString(myid);
            //PfClust::parse_domtabDir(tmpDir, newgroup, clean);
            break;
        }
        
        int big_chunk_index = this_chunk_info[0];
        int chunk_start = this_chunk_info[1];
        int chunk_end = this_chunk_info[2];
        int small_chunk_index = this_chunk_info[3];

        
        cout << "   " << myid << ": " << (chunk_end - chunk_start + 1) << " groups: "\
        << chunk_start << " to " << chunk_end << " from " << big_chunk_index << endl << std::flush;
        
        /* do the pipeline process for all the groups in this chunk */
        // input: fasta and group files for the big chunk
        string fastafile = inputDir + Utils::getPathSeparator() + Utils::intToString(big_chunk_index) + ".fasta";
        string groupfile = inputDir + Utils::getPathSeparator() + Utils::intToString(big_chunk_index) + ".group";
        
        /************************ read sequences and groups from input *********************************/
        
        map<string, string> seqMap, my_seqMap; // seqID => seq for all the sequences , big
        map<int, vector<string> > lenMap; // length => seqIDs with this length, big
        map<int, vector<string> > groupMap; // group => seqIDs, big
        vector<int> groupNames;
        
        //cout << "   " << myid << ": read fasta in memory" << endl << std::flush;
        /* read big fasta file */
        PfClust::read_fasta(fastafile, seqMap, lenMap);

        //cout << "   " << myid << ": read group in memory" << endl << std::flush;
        /* read group file */
        PfClust::read_group(groupfile, groupMap, groupNames);

        /************************ erase useless sequences to save memory *********************************/
        
        /* delete the groups that do not belong to this chunk */
        map<int, vector<string> >::iterator it;
        vector<string>::iterator seq_it;
        
        it = groupMap.find(groupNames.at(chunk_start));
        //cout << "   start -> " << it->first << endl;
        groupMap.erase(groupMap.begin(), it);
    
        it = groupMap.find(groupNames.at(chunk_end));
        if (it != groupMap.end()) {
            //cout << "   end -> " << it->first << endl;
            groupMap.erase(++it, groupMap.end());
        }
        
        /************************ print info for the groups in the node *********************************/
        
        // starting and ending group names
        cout << "   " << myid << ": " << seqMap.size() << " sequences, " << groupMap.size() << " groups: "\
        << groupMap.begin()->first << " to " << (--groupMap.end())->first << " from " << big_chunk_index << endl << std::flush;

        /*********************** generate domtab file for each group ***********************************/
        //cout << "   " << myid << ": start group sequences into fasta files" << endl << std::flush;
        //omp_set_num_threads(16);
        #pragma omp parallel for
        for (int chunk = chunk_start; chunk <= chunk_end; chunk++) {
            int gp_name = groupNames.at(chunk);
            int res = 0;
            PfClust::group_fasta_single(gp_name, groupMap.at(gp_name), tmpDir, tmpDir, seqMap); //generate fasta file for this chunk group
            //cout << "   " << myid << ": all fastas generated" << endl;
            
            string tmp_filename = tmpDir + Utils::getPathSeparator() + Utils::intToString(gp_name);
            string filename = outputDir + Utils::getPathSeparator() + Utils::intToString(gp_name);
            //cout << "   " << myid << ": current file: " << tmp_filename << endl;
            string msafile = tmp_filename + ".fasta";
            // if not singleton, align the sequences to create msa
            if ((groupMap.at(gp_name)).size() > 1) {
                res = 1;
                string mafft_cmd_prefix = mafft_cmd_path + " --auto --anysymbol --amino --quiet --thread 1";
                
                RunBioProg mafft(mafft_cmd_prefix, true, '>');
                mafft.setClean(clean);
                mafft.setInputfile(tmp_filename + ".fasta");
                mafft.setOutputfile(tmp_filename + ".mafft");
                res = mafft.runCmd();
                //cout << "   " << myid << " mafft command: " << mafft.getCmd() << endl << flush;
                msafile = tmp_filename + ".mafft";
                if (Utils::isFileEmpty(msafile)) {
                    cout << "   " << myid << " empty msa!! mafft command: " << mafft.getCmd() << endl;
                }
            }

            //cout << "   " << myid << ": all msas generated" << endl;
            // if mafft is run and run successfully, hmmbuild, build HMM for this group
            if ( res == 0 ) {
                string hmmbuild_cmd_prefix = hmmbuild_cmd_path + " --amino --cpu 1";
                
                // direct the output to /dev/null and name the hmm as the group name
                string options = "-o /dev/null -n " + Utils::intToString(gp_name);
                
                RunBioProg hmmbuild_cmd(hmmbuild_cmd_prefix, false, ' ');
                hmmbuild_cmd.setClean(clean);
                hmmbuild_cmd.setOption(options);
                hmmbuild_cmd.setInputfile(msafile);
                hmmbuild_cmd.setOutputfile(filename +  ".hmm");
                //cout << "   " << myid << " hmmbuild command: " << hmmbuild_cmd.getCmd() << endl;
                res = hmmbuild_cmd.runCmd();
                
                // if hmmbuild is successful
                if (res != 0 )
                    cout << "   " << myid << ": hmmbuild failed: " << hmmbuild_cmd.getCmd() << endl;
            }
            // print message about mafft command failing
            else
                cout << "   " << myid << ": mafft failed: " << endl;
            
        }
        /************************ parse domtab files *********************************/
        cout << "   " << myid << ": all hmms generated" << endl;
        
        cout << "   " << myid << ": done, chunk number " << small_chunk_index << endl;
        
        MPI_Send(0,0,MPI_INT,0,0,MPI_COMM_WORLD);
    }
}
