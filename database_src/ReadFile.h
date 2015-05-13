#ifndef READFILE_H
#define READFILE_H

/********************************************************************
 * ReadFile.h
 * Created by JJ Chai  on 05/13/2015.
 * Copyright (c) 2015 JJ Chai (ORNL). Allrights reserved.
********************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <map>
#include <unordered_map>
#include <sys/stat.h>
#include "utils.hpp"

typedef unsigned int UINT16;
using namespace std;

/* Read group file (seqID groupID seqLen) and save information */
void read_group_sg(const string & file_name, unordered_map<string, UINT16> & seq_gp, map<UINT16, vector<UINT16> > & group_len_map);

#endif //READFILE.HPP
