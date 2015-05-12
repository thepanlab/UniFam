#ifndef RUNBIOPROG_HPP
#define RUNBIOPROG_HPP

/********************************************************************
 * runBioProg.hpp
 
 * Class used to call a bioinformatics program using system call
 * Get the input and output file, with its own format, perform the system call
 * now only process 1 file
 
 * Created by JJ Chai on 10/25/2013.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
 ********************************************************************/

//should we do this dynamically???

#ifndef _WIN32
#include <unistd.h>
#endif


#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>

using namespace std;

class RunBioProg
{
public:
    /** Constructors.
     ** Create an instance by providing a command name (with full path if necessary).
     **/
    
    RunBioProg();
    RunBioProg(const string & );
    RunBioProg(const string & , const bool & , const char & );
    
    /** destructor **/
    virtual ~RunBioProg();
    
    /************** Options for the program **************/
    
    /* set the symbol for directing output */
    void setSymbol( const char & );
    
    /* set the order of input output */
    void setOrder ( const bool & );
    
    /* determine if the input file should be deleted afterwards */
    void setClean (const bool & );
    
    /* set the option for the command */
    void setOption( const string & );
    
    /************** Input and output for the program **************/
    void setInputfile( const string & );
    void setOutputfile( const string & );
    
    
    /************** full command retrieve and execution **************/
    /* get the full command */
    string getCmd();
    
    /* execute the command */
    int runCmd();
    
    
private:
    
    // name of the command (including path)
    string cmd_name;
    
    // order of input output, if input first then true
    bool input_first;
    
    // symbol for directing output
    char direct_symbol;
    
    // clean the input file afterwards?
    bool clean;
    
    
    // input and output files
    string inputfile, outputfile;
    
    //options
    string option;
    
};

#endif //RUNBIOPROG_HPP

