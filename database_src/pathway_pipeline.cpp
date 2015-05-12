#include <mpi.h> // for mpi
#include <vector>
#include "utils.hpp" // for useful utilities - namespace
#include <math.h>

/********************************************************
 * pathway_pipeline.cpp
 
 * pipeline to get annotation, and infer pathways using MetaCyc
 * to be parallelized
 
 * Created by JJ Chai  on 12/12/2013. Last modified 12/12/2013.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
 ********************************************************/

#define WORKTAG    1
#define DIETAG     2

using namespace std;

/* print usage */
void usage()
{
    cout << " [Usage]" << endl
    << "  pathway_pipeline [options] -cmd <python cmd file> -min <min_index> -max <max_index> " << endl
    << endl
    << " [Inputs]" << endl
    << "   -cmd <python cmd file> file with python pathway_pipeline command prefix" << endl
    << "   -ind <index file> file with indices of genomes to process" << endl
    << "   -min <min_index> starting index of genomes to process " << endl
    << "   -max <max_index> ending index of genomes to process " << endl
    << endl
    
    << "  [options]" << endl
    << "   -h/--help Display this help message" << endl
//    << "   -proc <procs> number of procs/threads per node, default: 1" << endl
    << endl;
}

///* Parse command line arguments */
///* pathway_pipeline [options] -cmd <python cmd file> -min <min_index> -max <max_index> */
//
void initializeArguments(int argc, char **argv,
                         string & python_cmd, // python pathway_pipeline command prefix
                         vector<int> & indices, // file include indices to process
                         int & min_index, // starting index of genomes to process
                         int & max_index) // ending index of genomes to process
{
    vector<string> Arguments;
    
    string pythonCmdFile = ""; // python command file
    string indexFile = ""; // index file
    python_cmd = "";
    min_index = -1;
    max_index = -1;
    while(argc--)
        Arguments.push_back(*argv++);
    
    for(int i = 1; i <= (int)Arguments.size()-1; i++)
    {
        // python command file
        if (Arguments[i] == "-cmd") {
            pythonCmdFile = Arguments[++i];
        }
        // index file
        else if (Arguments[i] == "-ind") {
            indexFile = Arguments[++i];
        }
        
        // starting index of genomes to process
        else if (Arguments[i] == "-min")
            min_index = Utils::stringToInt(Arguments[++i]);
        
        // ending index of genomes to process
        else if (Arguments[i] == "-max")
            max_index = Utils::stringToInt(Arguments[++i]);
        
        
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
    if ((pythonCmdFile== "") || (min_index == -1) || (max_index == -1))
    {
        cerr << "miss necessary parameter(s)"<<endl;
        cerr << "use -h/--help for usage" << endl;
        exit(1);
    }
    
    // read python command from file
    if(!Utils::isFileEmpty(pythonCmdFile)){
        ifstream mycmd(pythonCmdFile.c_str());
        getline(mycmd,python_cmd);
        mycmd.close();
    }
    
    // read genome indices to process
    if (!Utils::isFileEmpty(indexFile)) {
        ifstream myindex(indexFile.c_str());
        int ind;
        while (myindex >> ind) {
            indices.push_back(ind);
        }
        myindex.close();
    }
    
}

// master job
void MasterProcess(const int & min_index, const int & max_index)
{
    
    // total number of processes to be spawned
    int i, num_proc_spawn;
    
    int currentWorkID, num_proc, num_slave, result;
    
    int genomeNum = max_index - min_index + 1;
    
//    cout << "number of genomes is " << genomeNum << endl;
    
    // MPI calls and initiation
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
    
    num_slave = num_proc - 1; /* number of slave processes */
    
    num_proc_spawn = ((genomeNum <= num_slave) ? genomeNum : num_slave); /* number of processes to really spawn */
    
//    cout << "number of proces spawned: " << num_proc_spawn << endl;
    
    for (i = 0; i < num_proc_spawn; i++) {
        // index of the genomeSec working on now
        currentWorkID = i + min_index;
        // send the index to the slave process
        MPI_Send(&currentWorkID,1,MPI_INT, i+1, WORKTAG, MPI_COMM_WORLD);
        cout << "send to worker " << i+1 << " genome " << currentWorkID << endl;
    }
    
    if ( genomeNum > num_slave ) {
        // next file after each process gets a job
        currentWorkID = num_slave + min_index;
        
        while (currentWorkID <= max_index) {
            // receive message from any process that finishes its job
            MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            // send it a new job
            MPI_Send(&currentWorkID, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            cout << "send to worker " << status.MPI_SOURCE << " genome " << currentWorkID << endl;
            // another job is being done!
            currentWorkID++;
        }
    }
    
    // when all jobs are done, send signal everybody to quit (BROADCAST?)
    for (i = 1; i <= num_slave; i++) {
        MPI_Send(0, 0, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
    }
    //cout << "Master process is done." << endl;
}


//SlaveProcess(python_cmd);
// slave job
void SlaveProcess(const string & python_cmd)
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
        string sys_cmd = python_cmd + " --index " + Utils::intToString(currentWorkID);
        const int res = system(sys_cmd.c_str());
        if(res != 0)
            cout << "  *** Failed command " << sys_cmd << endl;
        MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    
}

// slave job
void SlaveProcess(const string & python_cmd, const vector<int> & indices)
{
    //cout << "use index file!" << endl;
    MPI_Status status;
    int currentWorkID, myid;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    // do jobs until master tells to stop
    while (true) {
        MPI_Recv(&currentWorkID, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        if (status.MPI_TAG == DIETAG) {
            break;
        }
        string sys_cmd = python_cmd + " --index " + Utils::intToString(indices.at(currentWorkID));
        const int res = system(sys_cmd.c_str());
        if(res != 0)
            cout << "  *** Failed command " << sys_cmd << endl;
        MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    
}
// main function

int main(int argc, char ** argv){
    
    string python_cmd;
    vector<int> indices;
    int min_index, max_index, myid;
    
    double start_time, finish_time, elapsed_time;
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    
    if (myid == 0) {
        start_time = MPI_Wtime();
        // display work start and time record
        cout << endl
        << "============================================================================"
        << endl << Utils::currentDateTime() << endl
        << " [Step 1] Initialize command line arguments -> " << endl << std::flush;
    }
    
    /* initialize command line arguments */
    initializeArguments(argc, argv, python_cmd, indices, min_index, max_index);
    int genomeNum = max_index - min_index + 1;
    
    
    if ( myid == 0 ) {
        
        cout << "python command prefix is: " << python_cmd << endl;
        cout << "process genome index from " << min_index << " to " << max_index << endl;
        cout << "Done!" << endl;
        cout << " [Step 2] Running python pipeline -> " << endl << std::flush;
        
        MasterProcess(min_index,max_index);
    }
    
    else
    {
        if (indices.size() == 0) {
            SlaveProcess(python_cmd);
        }
        else{
            SlaveProcess(python_cmd,indices);
        }
        
    }
    
    if( myid == 0 )
    {
        cout << " Done!" << endl;
        
        finish_time = MPI_Wtime();
        
        // display work end and time record
        cout << Utils::currentDateTime() << " All genomess are processed "<< endl
        << "Total Elapsed Time =  "
        << double(finish_time - start_time) << " [seconds]" << endl
        << "============================================================================"
        << std::endl << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}
