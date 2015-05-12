#include <mpi.h> // for mpi
#include <omp.h> // for omp
#include "utils.hpp" // for useful utilities - namespace
#include "directoryStructure.hpp" // for searching files in directory - class
#include "runBioProg.hpp" // for running system call of bioinfo programs - class
#include "pfclust.hpp"

/********************************************************
 * split_fasta.cpp
 
 * Split big fasta file into small chunks, preparation for running.
 * MPI version so that big data can be separated into many chunks.
 * 
 * Every big chunk will be further separated into smaller chunks depending
 * on the minimum sequence length in the big chunk. The information for each 
 * smaller chunk and the groups contained in it is written in the *.chunk 
 * and *.group files.
 *
 * This version generate fastas with sequences strictly from the groups, 
 * and does not consider the groups that these sequences can hit when hmmsearch
 * is used to find their new group again. For this purpose, use split_fasta_hmmsearch.cpp
 
 * Created by JJ Crosskey  on 05/12/2015. Last modified 05/12/2015.
 * Copyright (c) 2015 JJ Crosskey (ORNL). Allrights reserved.
 
 ********************************************************/

#define WORKTAG    1
#define DIETAG     2

using namespace std;

/* print usage */
/* pipeline [options] -f <fastafile> -g <groupfile> -n <number of chunks> -o <outputDir> -p <threads>*/
void usage()
{
    cout << " [Usage]" << endl
    << "  pipeline [options] -f <fastafile> -g <groupfile> -n <number of chunks> -o <outputDir> " << endl
    << endl
    << " [Inputs]" << endl
    << "   -f <fastafile> fasta file that includes all sequences" << endl
    << "   -g <groupfile> group file including per sequence group membership " << endl
    << "   -n <number of chunks> number of big chunks to split into" << endl
    << endl
    << " [Outputs]" << endl
    << "   -o <outputDir> directory for smaller files, and info for all chunks" << endl
    << endl
    << " [Options]" << endl
    << "   -h/--help Display this help message" << endl
    << "   -p <threads> number of threads every node, default: 1"  << endl
    << endl;
}

void initializeArguments(int argc, char **argv,
                         string & fastafile, // fasta file with all sequences to extract
                         string & groupfile, // group file
                         string & outputDir, // output directory for the group files
                         int & threads, // number of threads per node
                         int & num_big_chunks)// total number of big chunks to be splitted into
{
    vector<string> Arguments;
    
    fastafile  = "";
    groupfile = "";
    outputDir = "";
    num_big_chunks = 0;
    threads = 1;
    
    while(argc--)
        Arguments.push_back(*argv++);
    
    for(int i = 1; i <= (int)Arguments.size()-1; i++)
    {
        // input fastafile
        if (Arguments[i] == "-f") {
            fastafile = Arguments[++i];
        }
        // input groupfile
        else if(Arguments[i] == "-g")
            groupfile = Arguments[++i];
        
        // group file outputDir
        else if (Arguments[i] == "-o")
            outputDir = Arguments[++i];
        
        // number of big chunks
        else if (Arguments[i] == "-n")
            num_big_chunks = Utils::stringToInt(Arguments[++i]);
        
        // number of threads
        else if (Arguments[i] == "-p")
            threads = Utils::stringToInt(Arguments[++i]);
        
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
    if ((fastafile == "") || (groupfile == "") || (num_big_chunks == 0) || (outputDir == ""))
    {
        cerr << "miss necessary parameter(s)"<<endl;
        cerr << "use -h/--help for help" << endl;
        exit(1);
    }
    
}

// master job
void MasterProcess(const int & num_big_chunks, const int & threads)
{
    // total number of processes to be spawned
    int i, num_nodes_spawn, num_nodes, num_pieces;
    
    int currentWorkID, num_slave_nodes, result;

    cout << "number of big chunks is " << num_big_chunks << endl;
    
    // MPI calls and initiation
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD,&num_nodes);
    
    num_slave_nodes = num_nodes - 1; /* number of slave nodes */
    num_pieces = ceil((double)num_big_chunks/threads); // if every node works many threads, how many nodes are needed
    
    num_nodes_spawn = ((num_pieces <= num_slave_nodes) ? num_pieces : num_slave_nodes); /* number of processes to really spawn */
    
    /* every node gets several (num_threads) big chunks to work on */
    for (i = 1; i <= num_nodes_spawn; i++) {
        // index of the file working on now
        currentWorkID = i - 1;
        // send the index to the slave process
        MPI_Send(&currentWorkID, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
    }
    
    if ( num_pieces > num_slave_nodes ) {
        // next file after each process gets a job
        currentWorkID = num_slave_nodes ;
        
        while (currentWorkID < num_pieces) { /* currentWorkID starts at 0, num_big_chunks is the total number */
            // receive message from any process that finishes its job
            MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            // send it a new job
            MPI_Send(&currentWorkID, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            // another job is being done!
            currentWorkID++;
        }
    }
    
    // when all jobs are done, send signal everybody to quit (BROADCAST?)
    for (i = 1; i <= num_slave_nodes; i++) {
        MPI_Send(0, 0, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
    }
    //cout << "Master process is done." << endl;
}

// slave job
void SlaveProcess(const map< string, string > & seqMap, 
                  const map<int, vector<string> > & groupMap, const vector<int> & groupNames,
                  const int & gp_num, const int & num_groups, const int & threads, const int & num_big_chunks,
                  const string & outputDir)
{
    MPI_Status status;
    int currentWorkID, myid;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    // do jobs until master tells to stop
    while (true) {
        MPI_Recv(&currentWorkID, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        if (status.MPI_TAG == DIETAG) {
            break;
        }
	/* Each thread works on a big chunk, therefore each process works on multiple big chunks */
        int chunk_start = currentWorkID * threads;
        int chunk_end = min( (currentWorkID+1)*threads, num_big_chunks);
        
        cout << "   Process " << myid << ": working on chunks " << chunk_start << " to " << chunk_end << endl;
	int nProc = omp_get_max_threads();      /* maximum threads possible */
        omp_set_num_threads(threads > nProc ? nProc : threads); /* make sure thread count is within limit */
        
	/* Work on the big chunk with multiple threads */
        #pragma omp parallel for
        for (int j = currentWorkID * threads; j < chunk_end; j++) {
            
            int i, gp_start, gp_end, min_len, gp, chunks, chunk_group, len;
            vector<string> seqs;
            vector<string>::const_iterator it;
            
            gp_start = j * gp_num; //starting group number for this big chunk
            gp_end = min(num_groups-1, (j+1) * gp_num - 1); // ending group number (included in the chunk)
            
            
            /*************************** write subfiles *****************************/
            // file names
            string gp_file_base = outputDir + Utils::getPathSeparator() + Utils::intToString(j);
            string fastafilename = gp_file_base + ".fasta";
            string groupfilename = gp_file_base + ".group";
            
            if (Utils::isFileExist(fastafilename) && Utils::isFileExist(groupfilename) ) {
                cout << "   big chunk " << j << " already done" << endl;
            }
            else{
                /* find the maximum length and upper bound for this node */
                min_len = 5000;
                for (i = gp_start; i <= gp_end; i++) {
                    seqs = groupMap.at(groupNames.at(i));
                    for (it = seqs.begin(); it != seqs.end(); it++)
                    {
                        len = (int) seqMap.at(*it).length();
                        min_len = (len < min_len) ? len : min_len;
                    }
                    
                }
                
                ofstream myfasta(fastafilename.c_str()); /* fasta file with all sequences in this big chunk */
                ofstream mygroup(groupfilename.c_str()); /* group file with "seq group" format for all the sequences */
                for (gp = gp_start; gp <= gp_end; gp++) {
                    for (it = (groupMap.at(groupNames.at(gp))).begin(); it !=  (groupMap.at(groupNames.at(gp))).end(); it++) {
                        mygroup << *it << " " << groupNames.at(gp) << endl;
                        myfasta << ">" << *it << endl << seqMap.at(*it) << endl;
                    }
                }
                mygroup.close();
                myfasta.close();
            }
            
            /*************************** write chunk.info *****************************/
            string chunk_filename = gp_file_base + ".chunk";
            if (!Utils::isFileExist(chunk_filename)) {
                ofstream mychunk(chunk_filename.c_str());
                
                if (min_len/30 > 1) {
                    
                    int gps = gp_end - gp_start + 1; // actual number of groups in this big chunk
                    chunks = min(min_len/30,gps); // if min_len/30 is bigger than the number of groups, cannot have more than "gps" groups
                    if (gps % chunks == 0) {
                        chunk_group = gps/chunks;
                    }
                    else{
                        chunk_group = gps / (chunks - 1);
                        if (gps % (chunks - 1) == 0) {
                            chunks = chunks - 1;
                        }
                    }
                    
                    cout << "   big chunk " << j << " -> small chunks = " << chunks << endl;
                    for (i = 0; i < chunks; i++) {
                        mychunk << j << " " << i*chunk_group << " " << min(gps-1, (i+1)*chunk_group - 1) << endl;
                    }
                }
                // if shortest sequence is shorter than 60
                else{
                    cout << "   big chunk " << j << " -> small chunks = 1" << endl;
                    mychunk << j << " " << 0 << " " << (gp_end-gp_start) << endl;
                }
                
                mychunk.close();
            }
        }
        // signal master when done
        MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
        
    }
    
}


// main function
// function hmm_file prop count_file fastafile
int main(int argc, char ** argv){
    
    /********************** variable declaration *********************/
    
    string fastafile, groupfile, outputDir;
    int num_big_chunks, num_groups, gp_num, threads;
    map< string, string > seqMap;
    map<int, vector<string> > groupMap;
    vector<int> groupNames;
    
    /* initialize command line arguments */
    initializeArguments(argc, argv, fastafile, groupfile, outputDir, threads, num_big_chunks);
    Utils::mkdirIfNonExist(outputDir);
    
    int myid;
    double start_time, finish_time;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    
    /********************** display work start and time record *********************/
    if (myid == 0) {
        start_time = MPI_Wtime();
        cout << endl
        << "============================================================================"
        << endl << Utils::currentDateTime() << endl
        << " Beginning split_fasta" << endl
        << " [Step 1] Prepare splitting: Running -> " << endl << std::flush;
    }
    
    /********************** do stuff *********************/
    
    /* Read the sequences and groups into memory */
    PfClust::read_fastaseq(fastafile, seqMap);
    PfClust::read_group(groupfile, groupMap, groupNames);
    
    num_groups = (int) groupMap.size(); // total number of groups
    num_big_chunks = min(num_big_chunks,num_groups); // number of chunks cannot exceed number of groups
    /* number of groups every big chunk has */
    if (num_groups % num_big_chunks == 0) { // if every big chunk has same number of groups
        gp_num = num_groups/num_big_chunks;
    }
    else{ // otherwise, every chunk except the last one gets more groups than the last one
        gp_num = num_groups/(num_big_chunks - 1);
        // if in this case the last chunk gets nothing, decrease the number of big chunks
        if (num_groups % (num_big_chunks - 1) == 0) {
            num_big_chunks = num_big_chunks - 1;
        }
    }
    
    
    /* report numbers */
    if (myid == 0) {
        cout << "   " << num_groups << " total groups: "  << *(groupNames.begin()) << " to " << *(groupNames.rbegin()) << endl;
        cout << "   " << gp_num << " groups every big chunk" << endl;
        cout << "   " << num_big_chunks << " big chunks" << endl << endl;
        cout << Utils::currentDateTime() <<  " Done!" << endl;
        cout << " [Step 2] write files for all big chunks: Running -> " << std::flush;
        
        MasterProcess(num_big_chunks, threads);
    }
    
    else
    {
        SlaveProcess(seqMap, groupMap, groupNames, gp_num, num_groups, threads, num_big_chunks, outputDir);
    }
    
    if( myid == 0 )
    {
        cout << " Done!" << endl;
        
        finish_time = MPI_Wtime();
        
        // display work end and time record
        cout << Utils::currentDateTime() << " Ending split_fasta "<< endl
        << "Total Elapsed Time =  "
        << double(finish_time - start_time) << " [seconds]" << endl
        << "============================================================================"
        << std::endl << std::endl;
    }
    
    MPI_Finalize();

    return 0;
}
