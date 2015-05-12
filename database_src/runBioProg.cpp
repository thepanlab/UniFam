/********************************************************************
 * runBioProg.cpp
 
 * Class used to call a bioinformatics program using system call
 * Get the input and output file, with its own format, perform the system call
 * now only process 1 file
 
 * Created by JJ Chai on 10/25/2013. Last modified on 11/13/2013.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
 ********************************************************************/

#include "runBioProg.hpp"
#include "utils.hpp"

RunBioProg::RunBioProg()
{
	cmd_name = "";
    input_first = true;
    direct_symbol = ' ';
    clean = false;
    inputfile = "";
    outputfile = "";
    option = "";
    
}

RunBioProg::RunBioProg(const string & name)
{
	cmd_name = name;
    input_first = true;
    direct_symbol = ' ';
    clean = false;
    inputfile = "";
    outputfile = "";
    option = "";
    
}

RunBioProg::RunBioProg(const string & name, const bool & order, const char & symbol)
{
	cmd_name = name;
    input_first = order;
    direct_symbol = symbol;
    clean = false;
    inputfile = "";
    outputfile = "";
    option = "";
    
}

RunBioProg::~RunBioProg()
{

}

void RunBioProg::setSymbol(const char & symbol)
{
    direct_symbol = symbol;
}

void RunBioProg::setOrder(const bool &  order)
{
    input_first = order;
}

void RunBioProg::setClean(const bool & cleanup)
{
    clean = cleanup;
}

void RunBioProg::setOption(const string & Option)
{
    option = Option;
}

void RunBioProg::setInputfile(const string & input)
{
    inputfile = input;
}


void RunBioProg::setOutputfile(const string & output)
{
    outputfile = output;
}

string RunBioProg::getCmd()
{
    string full_cmd;
    full_cmd = cmd_name + " " + option;
    
    if (input_first)
    {
        full_cmd = full_cmd + " " + inputfile + " " + direct_symbol + " " + outputfile;
    }
    else
    {
        full_cmd = full_cmd + " " + outputfile + " " + direct_symbol + " " + inputfile;
    }
    
    return full_cmd;
}

/* run command, if output file already exists, overwrite it */
int RunBioProg::runCmd()
{
    string full_cmd = getCmd();
    const int res = system(full_cmd.c_str());
    if (res != 0) {
        cout << "   *** Failed command " << getCmd() << endl;
    }
    
    if (clean) {
        vector<string> inputfiles = Utils::commaStringToVector(inputfile, ' ');
        for (vector<string>::iterator it = inputfiles.begin(); it != inputfiles.end(); it++) {
            remove((*it).c_str());
        }
    }
    return res;
}
