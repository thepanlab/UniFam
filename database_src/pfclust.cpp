/********************************************************************
 * pfclust.cpp
 *
 * pfclust.cpp includes most functions in the project.
 *
 * Created by JJ Chai on 10/17/2013. Last modified 10/28/2013.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
 ********************************************************************/

#include "pfclust.hpp"

// characters to skip, until the EOF
#define MAX_WIDTH 10000

// maximum possible sequence length
#define MAX_SEQ_LENGTH 100000

namespace PfClust
{

    
    /********************************************************************/
    /********************************************************************/
//    // function to find the indices of the lower bound and upper bound
//    // lenMap is a map < length => seqIDs with this length>, keys are ordered by value
//    // does not need this function, because it's so simple....
//    void find_index(map<int, vector<string> > & lenMap, int & ub, int & lb)
//    {
//        map<int, vector<string> >::iterator itlow, itup;
//        
//        // 'up' points to the first element that is less than or equal to ub
//        // 'low' points to the first element that is less than (not equal) to lb
//        // if can't find, then they point to the end of the vector
//        
//        itlow = lenMap.lower_bound(lb);
//        itup = lenMap.upper_bound(ub);
//
//        
//    }
    
    /********************************************************************/
    /********************************************************************/
    // given a sequence length and the difference proportion, find the lower bound and upper bound
    void find_range(const int & length, const float & prop, int & ub, int & lb)
    {
        // if prop = 0, include sequences of all length,
        if (prop == 0) {
            ub = MAX_SEQ_LENGTH;
            lb = 0;
        }
        else
        {
            ub = ceil((float) length / prop);
            lb = floor((float) length * prop);
        }
    }
    
    /********************************************************************/
    /********************************************************************/
    
    // find the length of hmm from .hmm file
    // if the file is empty, return 0 as the length
    // if the file is not correct, the result is not a number, also return 0
    // 0 should be a red flag! (something is wrong)
    
    int get_hmm_len(const string & hmmfile)
    {
        string line;
        int length;
        
        // read the sequence length file if it exists and not empty
        if (!Utils::isFileEmpty(hmmfile)) {
            
            ifstream myfile;
            myfile.open(hmmfile.c_str());
            
            // length of hmm is on the 3rd line
            for (int lineno = 0 ; lineno < 3; lineno++)
            {
                if (lineno==2)
                {
                    myfile >> line;
                    myfile >> length; // get the length (number of nodes) of this hmm
                }
                else
                    getline(myfile,line);
            }
            myfile.close();
            
            return length;
        }
        else
        {
            std::cerr << "hmmfile " << hmmfile << " is empty";
            return 0;
        }
    }
    
    /********************************************************************/
    /********************************************************************/
    // get all seqIDs from group_sg, and fine the number of sequences with bit score higher than the cutoff
    void get_seqIDs(const string & group_sg, map<string, int> & seqID_index, const float & cutoff, int & num_seqs) // group file, per sequence
    {
        // check if file is empty or does not exist
        if (Utils::isFileEmpty(group_sg)) {
            Utils::exitWithError("group file is empty \n");
        }
        else{
            string seqID, groupID;
            float bitscore;
            int i = 0;
            num_seqs = 0;
            ifstream mygroup(group_sg.c_str());
            
            while (mygroup >> seqID) {
                mygroup >> groupID >> bitscore;
                i++;
                
                if (bitscore >= cutoff) {
                    num_seqs++;
                    seqID_index[seqID] = num_seqs;
                }
                mygroup.ignore(MAX_WIDTH,'\n');
            }
            
            mygroup.close();
            
            cout << "Out of " << i << " " << num_seqs << "( " << seqID_index.size()
            << ") has bit score greater than " << cutoff << endl;
        }
    }
    
    /********************************************************************/
    /********************************************************************/
    
    // read fasta file and store sequences in a map <seqID => seq>,
    // lenMap : <length => sequence IDs with sequences of this length>
    void read_fasta(const string & fastafile, map<string, string> & seqMap, map<int, vector<string> > & lenMap)
    {
        string line, seqID, seq;
        int seqLen;
        char first_char;
        seq = "";
        
        // check if file is empty or does not exist
        if (Utils::isFileEmpty(fastafile)) {
            Utils::exitWithError("fasta file is empty \n");
        }
        
        else
        {
            //read the fasta file
            ifstream myfile(fastafile.c_str());
            
            while (myfile.get(first_char)) {
                
                // if the line contains sequence name
                if ( first_char == '>')
                {

                    
                    // if this is not the first sequence in the file, build map
                    if (seq.compare("")!=0) {
                        seqMap[seqID] = seq;
                        seqLen = (int) seq.length();
                        lenMap[seqLen].push_back(seqID);
                        // get sequence ID until the first space, and ignore the rest of annotation
                        myfile >> seqID;
                        myfile.ignore(MAX_WIDTH,'\n');
                    }
                    else{
                        // get sequence ID until the first space, and ignore the rest of annotation
                        myfile >> seqID;
                        myfile.ignore(MAX_WIDTH,'\n');
                    }
                    
                    //reset seq
                    seq = "";
                }
                else
                {
                    myfile.unget();
                    getline(myfile,line);
                    seq = seq + line;
                }
            }
            
            // the last sequence
            seqMap[seqID] = seq;
            seqLen = (int) seq.length();
            lenMap[seqLen].push_back(seqID);
            
            myfile.close();
        }
        
    }
    
    
    //read fasta file and store sequences in a map <seqID => seq>
    void read_fastaseq(const string & fastafile, map<string, string> & seqMap)
    {
        string line, seqID, seq;
        char first_char;
        seq = "";
        
        // check if file is empty or does not exist
        if (Utils::isFileEmpty(fastafile)) {
            Utils::exitWithError("fasta file is empty \n");
        }
        
        else
        {
            //read the fasta file
            ifstream myfile(fastafile.c_str());
            
            while (myfile.get(first_char)) {
                
                // if the line contains sequence name
                if ( first_char == '>')
                {
                    // if this is not the first sequence in the file, build map
                    if (seq.compare("")!=0) {
                        seqMap[seqID] = seq;
                        // get sequence ID until the first space, and ignore the rest of annotation
                        myfile >> seqID;
                        myfile.ignore(MAX_WIDTH,'\n');
                    }
                    else{
                        // get sequence ID until the first space, and ignore the rest of annotation
                        myfile >> seqID;
                        myfile.ignore(MAX_WIDTH,'\n');
                    }
                    
                    //reset seq
                    seq = "";
                }
                else
                {
                    myfile.unget();
                    getline(myfile,line);
                    seq = seq + line;
                }
            }
            
            // the last sequence
            seqMap[seqID] = seq;
            myfile.close();
        }
        
    }
    
    /********************************************************************/
    /********************************************************************/
    
    // read a (big) fasta file and write their length information to a file
    void get_fasta_len(const string & fastafile, const string &lengthfile)
    {
        string line, seqID, seq;
        map<string, int> seqLenMap;
        char first_char;
        seq = "";
        
        // check if file is empty or does not exist
        if (Utils::isFileEmpty(fastafile)) {
            Utils::exitWithError("fasta file is empty \n");
        }
        
        else
        {
            //read the fasta file
            ifstream myfile(fastafile.c_str());
            
            while (myfile.get(first_char)) {
                
                // if the line contains sequence name
                if ( first_char == '>')
                {
                    
                    // if this is not the first sequence in the file, build map
                    if (seq.compare("")!=0) {
                        seqLenMap[seqID] = (int) seq.length();
                        // get sequence ID until the first space, and ignore the rest of annotation
                        myfile >> seqID;
                        myfile.ignore(MAX_WIDTH,'\n');
                    }
                    else{
                        // get sequence ID until the first space, and ignore the rest of annotation
                        myfile >> seqID;
                        myfile.ignore(MAX_WIDTH,'\n');
                    }
                    
                    //reset seq
                    seq = "";
                }
                else
                {
                    myfile.unget();
                    getline(myfile,line);
                    seq = seq + line;
                }
            }
            
            // the last sequence
            seqLenMap[seqID] = (int) seq.length();
            
            myfile.close();
            
            // write sequence lengths to file
            ofstream lengthout(lengthfile.c_str());
            for (map<string, int>::iterator it = seqLenMap.begin(); it != seqLenMap.end(); it++) {
                lengthout << it->first << "\t" << it->second << endl;
            }
            
            lengthout.close();
            
        }

    }
    
    // read a (big) fasta file and store their length information in a vector
    void read_fasta_len(const string & fastafile, map<string, int> & seqLenMap)
    {
        string line, seqID, seq;
        char first_char;
        seq = "";
        
        // check if file is empty or does not exist
        if (Utils::isFileEmpty(fastafile)) {
            Utils::exitWithError("fasta file is empty \n");
        }
        
        else
        {
            //read the fasta file
            ifstream myfile(fastafile.c_str());
            
            while (myfile.get(first_char)) {

                // if the line contains sequence name
                if ( first_char == '>')
                {
                    
                    // if this is not the first sequence in the file, build map
                    if (seq.compare("")!=0) {
                        seqLenMap[seqID] = (int) seq.length();
                        // get sequence ID until the first space, and ignore the rest of annotation
                        myfile >> seqID;
                        myfile.ignore(MAX_WIDTH,'\n');
                    }
                    
                    else{
                        // get sequence ID until the first space, and ignore the rest of annotation
                        myfile >> seqID;
                        myfile.ignore(MAX_WIDTH,'\n');
                    }
                    
                    //reset seq
                    seq = "";
                }
                else
                {
                    // get sequence ID until the first space, and ignore the rest of annotation

                    myfile.unget();
                    getline(myfile,line);
                    seq = seq + line;
                }
            }
            
            // the last sequence
            seqLenMap[seqID] = (int) seq.length();
            
            myfile.close();
        }
    }
    /********************************************************************/
    /********************************************************************/
    //given the filename of hmm, do hmmsearch
    
    void do_hmmsearch(const string & hmmFilename, // HMM file name
                      const string & outputDir, // output directory
                      map<string, string> & seqMap, // seqID => seq map
                      const map<int, vector<string> > & lenMap, // length => seqIDs with this length
                      const float & prop, // length cutoff for shorter seq(hmm) / longer seq(hmm)
                      const string & hmmsearch_cmd_prefix, // first part of hmmsearch command
                      const bool & log, // save the log file?
                      const bool & clean) // clean up HMM file after the search?
    {
        
        int ub,lb,length; // upper bound and lower bound for sequence lengths to be included in the fasta file
        string filename, options;
        
        /* if outputDir is not a directory, make a directory with that name */
        Utils::mkdirIfNonExist(outputDir);
        
        filename = Utils::getFilename2(hmmFilename);
        
        /* domtab file name -- output */
        string dom_tab_file = outputDir + Utils::getPathSeparator() + filename + ".domtab";
        
        /************** check if the search is already done ***********/
        
        bool search_done = false;
        
//        if(Utils::isFileExist(dom_tab_file)){
//            
//            string lastline = Utils::last_line(dom_tab_file);
//            
//            // existing hmmsearch was not done correctly
//            if (lastline.compare("# [ok]")==0)
//                search_done = true;
//
//        }
        
        /************* search if job is not already done before *************/

        if(!search_done){
            
            length = get_hmm_len(hmmFilename);
            
            /* check if hmm file is good */
            if (length == 0)
            {
                cout << " hmm file " << hmmFilename << " is empty or does not exist, or does not have right format" << endl;
            }
            else
            {
                /* construct hmmsearch command */
                RunBioProg my_hmmsearch(hmmsearch_cmd_prefix, true, ' ');
                my_hmmsearch.setClean(clean);
                
                string fasta_file = outputDir + Utils::getPathSeparator() + filename + ".fa";
                /****** find candidate sequences *******/
                
                // find the lower bound and upper bound from the sequence length and threshold
                find_range(length,prop,ub,lb);
                
                map<int, vector<string> >::const_iterator itlow,itup, itfor;
                vector<string>::iterator itstr;

                itlow = lenMap.lower_bound(lb);
                itup = lenMap.upper_bound(ub);
                
                /**** save the candidate sequences in a fasta file *****/
                
                ofstream myfasta;
                myfasta.open(fasta_file.c_str());
                
                if(myfasta.fail()){
                    cout << "fasta file cannot be opened to write\n";
                }
                
                // write all sequences with lengths within range to fasta file
                for (itfor = itlow; itfor!=itup; itfor++) {
                    vector<string> seqIDs = itfor->second;
                    //for each sequence with this length:
                    for (itstr = seqIDs.begin(); itstr != seqIDs.end(); itstr++) {
                        myfasta << ">" << *itstr << endl
                        << seqMap[*itstr] << endl;
                    }
                }
                
                myfasta.close();
                /************* do hmmsearch ****************/
                //cout << "do hmmsearch for " << filename << ".hmm" << endl;
                
                options = "--domtblout " +  dom_tab_file;
                if (log)
                    options = "-o " + outputDir + Utils::getPathSeparator() + filename + ".log " + options;
                else
                    options = "-o /dev/null " +  options;
                
                my_hmmsearch.setOption(options);
                // both files are "input", no output specified except --domtab domtab_file
                my_hmmsearch.setInputfile(hmmFilename + " " + fasta_file);
                
                //cout << my_hmmsearch.getCmd() << endl;
                my_hmmsearch.runCmd();
                
                /******** remove the fasta files *********/
                remove(fasta_file.c_str());
            }
        }
        
        /**** search is already done correctly ****/
        else
        {
            cout << dom_tab_file << " already exists!\n";
        }
        
    }
    
    // parse the group file (with bit score attached)
    void parse_group(const string & filename, // group file name, with bit score
                     map<string, group_score> & scores, // resulting map
                     const bool & clean) // remove the group file after parsing??
    {
        if (Utils::isFileEmpty(filename)){
            cout << "group file "<< filename << " is empty \n";
        }
        else{
            string seqID;
            group_score group_bitscore;
            
            ifstream mygroup(filename.c_str());
            
            while (mygroup >> seqID) {
                mygroup >> group_bitscore.groupID >> group_bitscore.bitscore;
                mygroup.ignore(MAX_WIDTH,'\n');
                scores[seqID] = group_bitscore;
            }
            mygroup.close();
            
            if(clean)
                remove(filename.c_str());
        }
    }
    
    //parse all the groups in the directory, and generate overall group file
    void parse_groupDir(const string & groupfileDir,
                        const string & output,
                        const bool & clean)
    {
        // stores all files in the directory
        vector<string> groupFilenames;
        // get all domtab files in the input directory
        Utils:: getFiles(groupfileDir, "*", groupFilenames);

        map<string, group_score> scores;
        map<string, group_score>::iterator map_it;
        string filename, seqID;
        float best_score;
        group_score group_bitscore;
        vector<string>::const_iterator it = groupFilenames.begin();
        
        // parse the first groupfile
        parse_group(*it, scores, clean);
        
        for (it = ++(groupFilenames.begin()); it != groupFilenames.end(); it++) {
            filename = *it;
            if (Utils::isFileEmpty(filename)){
                cout << "group file "<< filename << " is empty \n";
            }
            else{
                ifstream mygroup(filename.c_str());
                
                while (mygroup >> seqID) {
                    mygroup >> group_bitscore.groupID >> group_bitscore.bitscore;
                    mygroup.ignore(MAX_WIDTH,'\n');
                    
                    // see if the sequence is already in the map
                    map_it = scores.find(seqID);
                    if (map_it == scores.end()) {
                        scores[seqID] = group_bitscore;
                    }
                    // if it is, compare the current best score with the new one
                    else{
                        best_score = (scores.at(seqID)).bitscore;
                        if (group_bitscore.bitscore > best_score)
                            scores[seqID] = group_bitscore;
                    }
                }
                mygroup.close();
                
                if(clean)
                    remove(filename.c_str());
            }
        }
        
        // print overall group membership in the new file
        ofstream myout(output.c_str());
        for (map_it = scores.begin(); map_it != scores.end(); map_it++) {
            myout << map_it->first << "\t" << (map_it->second).groupID << "\t" << (map_it->second).bitscore << endl;
        }
        myout.close();
    }
    /********************************************************************/
    /********************************************************************/
    
    //parse all the groups in the directory, and generate set covering input file
    void setGen_groupDir(const string & groupfileDir, const string & output,
                         const string & groupIDfile,
                         const map<string, int> & seqID_index,
                         const int & num_seqs,
                         const float & cutoff,
                         const bool & clean)
    {
        // stores all files in the directory
        vector<string> groupFilenames;
        // get all domtab files in the input directory
        Utils:: getFiles(groupfileDir, "*", groupFilenames);
        
        map<int, vector<int> > sets;
        map<int, vector<int> >::iterator set_it;
        vector<int>::iterator seq_it;
        string filename, seqID;
        int groupID;
        float bitscore;
        
        vector<string>::const_iterator it;
        
        for (it = groupFilenames.begin(); it != groupFilenames.end(); it++) {
            filename = *it;
            if (Utils::isFileEmpty(filename)){
                cout << endl << "group file "<< filename << " is empty \n";
            }
            else{
                ifstream mygroup(filename.c_str());
                
                while (mygroup >> seqID) {
                    mygroup >> groupID >> bitscore;
                    mygroup.ignore(MAX_WIDTH,'\n');
                    
                    // if the bitscore is bigger than the cutoff, put it in the set
                    if (bitscore >= cutoff) {
                        sets[groupID].push_back(seqID_index.at(seqID));
                    }
                }
                mygroup.close();
                
                if(clean)
                    remove(filename.c_str());
            }
        }
        
        cout << sets.size() << " groups are included in the file. cutoff = " << cutoff << endl;
        // print set information file
        ofstream myout(output.c_str());
        ofstream groupout(groupIDfile.c_str());
        myout << num_seqs << " " << sets.size() << endl;
        // for each set
        for (set_it = sets.begin(); set_it != sets.end(); set_it++) {
            if ((set_it->second).size() > 0) {
                groupout << set_it->first << endl;
                
                // print all genes in the set
                for (seq_it = (set_it->second).begin(); seq_it != (set_it->second).end(); seq_it++) {
                    myout << *seq_it << " ";
                }
                myout << endl;
            }

        }
        myout.close();
        groupout.close();
    }
    
    /********************************************************************/
    /********************************************************************/
    // parse the domain wise table from hmmsearch, find the indices and bit score of the target sequences
    // build map < seqID => bitwise score per residue >
    void parse_domtab(const string & filename, // domtab file name
                          map<string,float> & targetScores,
                      const bool & clean) //resulting map
    {
        // if file is empty, print message, but don't quit
        if (Utils::isFileEmpty(filename))
        {
            cout << "domain tab file "<< filename << " is empty \n";
        }
        
        //if not either, parse the file
        else
        {
            //myfile -- infile stream
            ifstream myfile(filename.c_str());
            
            char first_char; // first character of the line
            string line; //line
            string lastSeq = " "; // previous target sequence ID
            
            //vector<string> targetSeqs;
            
            while(myfile.get(first_char))
            {
                // ignore lines that start with "#"
                if (first_char == '#') {
                    myfile.ignore(MAX_WIDTH,'\n');
                }
                
                // target sequences
                else
                {
                    myfile.unget(); //put the first character back in the stream
                    
                    string targetSeqID, targetLen, targetScore, others;
                    
                    myfile >> targetSeqID; //get seqID
                    
                    // if it's the same sequence's another domain, ignore; otherwise, update vectors
                    
                    if (lastSeq.compare(targetSeqID)!=0)
                    {
                        myfile >> others >> targetLen >> others >> others >> others >> others >> targetScore;
                        
                        // ignore the rest of the line after the ID, length, and total score are read
                        myfile.ignore(MAX_WIDTH,'\n');
                        
                        // put the new sequence and its bitwise score per residue in the vector
                        if ((targetScore.compare("")) != 0 && (targetLen.compare("")!=0) ) {
                            
                            targetScores[targetSeqID] = Utils::stringToFloat(targetScore) / Utils::stringToFloat(targetLen);
                            
                            //update the last target sequence
                            lastSeq = targetSeqID;
                            
                        }
                        else
                        {
                            // print warning, and do not save in the vector
                            cout << "*** Warning: in file " << filename << "there is something wrong with format, please check... \n";
                        }
                        
                        
                        //check, for debugging
/*                     cout << "targetSeqID is " << targetSeqID << "\t" << "length is " << targetLen << "\t" << "score is " 
 * 			    << targetScores[targetSeqID] << endl;
 */
                    }
                    else
                        myfile.ignore(MAX_WIDTH,'\n');
                }
            }// end reading file
            
            myfile.close();
            if (clean) {
                remove(filename.c_str());
            }
        }
    }
    
    typedef map<string,float> hmmMap;
    typedef map<string, hmmMap> seqScoreMap;

    // parse all domtab files in a directory
    void parse_domtabDir(const string & domtabDir, const string & groupfile, const bool & clean)
    {
        // stores all files in the directory
        vector<string> domtabFilenames;
        // get all domtab files in the input directory
        Utils:: getFiles(domtabDir, ".domtab", domtabFilenames);
        
        // seqID => < hmm -> score>
        seqScoreMap bitscore_allseqs;
        
        // parse domain tab files one by one
        // i: hmm index
        for (vector<string>::iterator it = domtabFilenames.begin(); it != domtabFilenames.end(); it++) {
            // hmm -> score map
            hmmMap targetScores;
            
            //string domtabfilename = *it;
            // index/name of the domtab file
            string domtabIndex = Utils::getFilename2(*it);
            
            // store the target sequence indices and the scores
            parse_domtab(*it, targetScores, clean);
            
            // if the domtab includes some target sequences
            
            if (targetScores.size() != 0) {
                
                //cout << *it << "\n";
                
                for (hmmMap::iterator hmmit = targetScores.begin(); hmmit != targetScores.end(); hmmit++) {
                    (bitscore_allseqs[hmmit->first])[domtabIndex] = hmmit->second;
                }
            }
            else
            {
                cout << *it << ": \t nothing found! \n **********";
            }
            
        }
        
        ofstream mygroup(groupfile.c_str());
        
        // find the best matching hmm for each sequence, max index and value from array
        
        for (seqScoreMap::iterator seqit=bitscore_allseqs.begin(); seqit != bitscore_allseqs.end(); seqit++) {
            
            hmmMap::iterator it;
            
            float max_val = (seqit->second).begin()->second;
            string max_hmm = (seqit->second).begin()->first;
            
            for (it = (seqit->second).begin(); it != (seqit->second).end(); it++) {
                if (it->second > max_val) {
                    max_val = it->second;
                    max_hmm = it->first;
                }
            }
            
            /** print result for this sequence in the group file **/
            
            //left adjust, and fixed width
            mygroup.flags(std::ios::left);
            
            mygroup.width(30);
            mygroup << seqit->first;
            
            mygroup.width(20);
            mygroup <<  max_hmm; //offset of model number
            
            mygroup.width(20);
            mygroup <<  max_val << endl;
        }
        mygroup.close();
    }
    
    
    // parse all domtab files in a directory
    void parse_domtabs(const vector<string> & domtabFilenames, const string & groupfile, const bool & clean)
    {
        
        // seqID => < hmm -> score>
        seqScoreMap bitscore_allseqs;
        
        // parse domain tab files one by one
        // i: hmm index
        for (vector<string>::const_iterator it = domtabFilenames.begin(); it != domtabFilenames.end(); it++) {
            // hmm -> score map
            hmmMap targetScores;
            
            //string domtabfilename = *it;
            // index/name of the domtab file
            string domtabIndex = Utils::getFilename2(*it);
            
            // store the target sequence indices and the scores
            parse_domtab(*it, targetScores, clean);
            
            // if the domtab includes some target sequences
            
            if (targetScores.size() != 0) {
                
                //cout << *it << "\n";
                
                for (hmmMap::iterator hmmit = targetScores.begin(); hmmit != targetScores.end(); hmmit++) {
                    (bitscore_allseqs[hmmit->first])[domtabIndex] = hmmit->second;
                }
            }
            else
            {
                cout << *it << ": \t nothing found! \n **********";
            }
            
        }
        
        ofstream mygroup(groupfile.c_str());
        
        // find the best matching hmm for each sequence, max index and value from array
        
        for (seqScoreMap::iterator seqit=bitscore_allseqs.begin(); seqit != bitscore_allseqs.end(); seqit++) {
            
            hmmMap::iterator it;
            
            float max_val = (seqit->second).begin()->second;
            string max_hmm = (seqit->second).begin()->first;
            
            for (it = (seqit->second).begin(); it != (seqit->second).end(); it++) {
                if (it->second > max_val) {
                    max_val = it->second;
                    max_hmm = it->first;
                }
            }
            
            /** print result for this sequence in the group file **/
            
            //left adjust, and fixed width
            mygroup.flags(std::ios::left);
            
            mygroup.width(30);
            mygroup << seqit->first;
            
            mygroup.width(20);
            mygroup <<  max_hmm; //offset of model number
            
            mygroup.width(20);
            mygroup <<  max_val << endl;
        }
        mygroup.close();
    }
    
    // parse domtab file generated from hmmsearch, with multiple hmms in the hmm database, instead of only one
    void parse_domtabfile(const string & domtabfile, const string & groupfile)
    {
        
        // seqID => < hmm -> score>
        seqScoreMap bitscore_allseqs;
        
        ifstream myfile(domtabfile.c_str());
        
        char first_char; // first character of the line
        string line; //line
        string lastSeq = " "; // previous target sequence ID
        
        //vector<string> targetSeqs;
        
        while(myfile.get(first_char))
        {
            // ignore lines that start with "#"
            if (first_char == '#') {
                myfile.ignore(MAX_WIDTH,'\n');
            }
            
            // target sequences
            else
            {
                myfile.unget(); //put the first character back in the stream
                
                string targetSeqID, targetLen, targetScore, others, groupID;
                int domain;
                
                myfile >> targetSeqID >> others >> targetLen >> groupID >> others >> others >> others >> targetScore >> others >> domain; //get seqID
                // ignore the rest of the line after the ID, length, and total score are read
                myfile.ignore(MAX_WIDTH,'\n');
                
                // if it's the same sequence's another domain, ignore; otherwise, update vectors
                if (domain == 1)
                {
                    
                    // put the new sequence and its bitwise score per residue in the vector
                    if ((targetScore.compare("")) != 0 && (targetLen.compare("")!=0) ) {
                        
                        bitscore_allseqs[targetSeqID][groupID] = Utils::stringToFloat(targetScore) / Utils::stringToFloat(targetLen);
                    }
                    else
                    {
                        // print warning, and do not save in the vector
                        cout << "*** Warning: in file " << domtabfile << "there is something wrong with format, please check... \n";
                    }
                    
                    
                    //check, for debugging
/*                     cout << "targetSeqID is " << targetSeqID << "\t" << "length is " << targetLen << "\t" << "score is " \
 *                     <<  bitscore_allseqs[targetSeqID][groupID] << endl;
 */
                }
            }
        }// end reading file

        
        myfile.close();

        
        ofstream mygroup(groupfile.c_str());
        
        // find the best matching hmm for each sequence, max index and value from array
        
        for (seqScoreMap::iterator seqit=bitscore_allseqs.begin(); seqit != bitscore_allseqs.end(); seqit++) {
            
            hmmMap::iterator it;
            
            float max_val = (seqit->second).begin()->second;
            string max_hmm = (seqit->second).begin()->first;
            
            for (it = (seqit->second).begin(); it != (seqit->second).end(); it++) {
                if (it->second > max_val) {
                    max_val = it->second;
                    max_hmm = it->first;
                }
            }
            
            /** print result for this sequence in the group file **/
            
            //left adjust, and fixed width
            mygroup.flags(std::ios::left);
            
            mygroup.width(30);
            mygroup << seqit->first << "\t";
            
            mygroup.width(20);
            mygroup <<  max_hmm << "\t"; //offset of model number
            
            mygroup.width(20);
            mygroup <<  max_val << endl;
        }
        mygroup.close();
    }
    
    // reorder group file by sequence ID
    void reorder_group(const string & filename, string & ordered_file)
    {
        // check if input file exists and is not empty
        if(Utils::isFileEmpty(filename))
        {
            Utils::exitWithError("Group file is empty.\n");
        }
        else
        {
            // map sequence to group
            map<string, string> groups;
            string seqID, group;
            
            ifstream infile(filename.c_str());
            
            
            while (infile >> seqID) {
                infile >> group;
                infile.ignore(MAX_WIDTH,'\n');
                groups[seqID] = group;
            }
            infile.close();
        
            ofstream outfile(ordered_file.c_str());
            // output file sorted by sequence ID
            for (map<string,string>::iterator it = groups.begin(); it != groups.end(); it++) {
                
                outfile.flags(std::ios::left);
                outfile.width(30);
                outfile << it->first;
                outfile.width(20);
                outfile << it->second << endl;
            }
            outfile.close();
        }
    }
    
    // generate group file by group from group file by sequence
    void sg_to_gs(const string & group_sg, const string & group_gs)
    {
        if(Utils::isFileEmpty(group_sg))
        {
            Utils::exitWithError("Group file is empty.\n");
        }
        else
        {
            // group => <seq1, seq2, ... seqn>
            map< string, vector<string> > groups;
            string seqID, group;
            vector<string> groupVec;
            
            ifstream infile(group_sg.c_str());
            
            while (infile >> seqID) {
                //infile >> seqID >> group;
                infile >> group;
                //cout << seqID << "\t" << group << endl;
                infile.ignore(MAX_WIDTH,'\n');
                groups[group].push_back(seqID);
                //cout << "push sequence " << seqID << endl;

            }
            infile.close();
            
            // write output file
            ofstream outfile(group_gs.c_str());
            for (map< string,vector<string> >::iterator it = groups.begin(); it != groups.end(); it++)
            {
                //left adjust
                outfile.flags(std::ios::left);
                
                // group name
                outfile.width(20);
                outfile << it->first << ": : " ;
                
                groupVec = it->second;
                
                for (vector<string>::iterator seq_it = groupVec.begin(); seq_it!=groupVec.end(); seq_it++)
                {
                    outfile << *seq_it << " ";
                }
                
                outfile << " : " << groupVec.size() << endl;
            
            }
            
            outfile.close();
        }
    }
    
    // read group information into memory <groupID, seqIDs in the group>
    void read_group(const string & group_sg, map<int, vector<string> > & groupMap, vector<int> & groupNames)
    {
        if(Utils::isFileEmpty(group_sg))
        {
            Utils::exitWithError("Group file is empty.\n");
        }
        else
        {
            string seqID;
            int group;
            
            ifstream infile(group_sg.c_str());
            
            while (infile >> seqID) {
                //infile >> seqID >> group;
                infile >> group;
                //cout << seqID << "\t" << group << endl;
                infile.ignore(MAX_WIDTH,'\n');
                groupMap[group].push_back(seqID);
                //cout << "push sequence " << seqID << endl;
                
            }
            infile.close();
            
            /* save the group names in a vector */
            map<int, vector<string> >::iterator it;
            for (it = groupMap.begin(); it != groupMap.end(); it++) {
                groupNames.push_back(it->first);
            }
        }
    }
    
    
    // generate fasta files, one for each group
    void group_fasta(const map< int, vector<string> > & groups, // group map <group => seqs in the group>
                     const string & nonSingletonDir, // output directory for all the groups that are not singletons
                     const string & singletonDir, // output directory for all the singleton groups
                     map<string, string> & seqMap) // sequence map <seqID, seq>
    {

        if (!Utils::isDirectory(nonSingletonDir)) {
            
            string cmd = "mkdir -p "+ nonSingletonDir;
            const int res = system(cmd.c_str());
            if (res != 0 )
                Utils::exitWithError("*** Error: Failed command: " + cmd);
        }
        
        if (!Utils::isDirectory(singletonDir)) {
            string cmd = "mkdir -p "+ singletonDir;
            const int res = system(cmd.c_str());
            if (res != 0 )
                Utils::exitWithError("*** Error: Failed command: " + cmd);
        }
        
        /* for reach group, generate the files */
        
        // name of the fasta file
        string fastafile;
        
        for (map<int, vector<string> >::const_iterator it = groups.begin(); it != groups.end(); it++)
        {
            //cout << it->first << endl;
            //singleton group
            if ((it->second).size() == 1)
            {
                fastafile = singletonDir + "/" +  Utils::intToString(it->first)  + ".fasta";
            }
            
            //non-singleton group
            else
            {
                fastafile = nonSingletonDir + "/" + Utils::intToString(it->first) + ".fasta";
            }
            
            ofstream fout(fastafile.c_str());
            
            for (vector<string>::const_iterator seq_it = (it->second).begin(); seq_it != (it->second).end(); seq_it++)
            {
                //cout << ">" << *seq_it << endl;
                fout << ">" << *seq_it << endl; /* write the seqID */
                
                fout << seqMap[*seq_it] << endl; /* write the actual sequence */
            }
            
            fout.close();
            
        }
        
    }
    // generate fasta file for a single group
    void group_fasta_single(const int & groupID, const vector<string> & seqIDs, // group map <group => seqs in the group>
                     const string & nonSingletonDir, // output directory for all the groups that are not singletons
                     const string & singletonDir, // output directory for all the singleton groups
                     map<string, string> & seqMap) // sequence map <seqID, seq>
    {
        string fastafile;
        
        // if group has no sequence, throw warning
        if (seqIDs.size() < 1) {
            cout << "*** warning: " << groupID << " does not have sequences! " << endl;
        }
        
        else
        {
            //singleton group
            if (seqIDs.size() == 1)
            {
                fastafile = singletonDir + Utils::getPathSeparator() +  Utils::intToString(groupID)  + ".fasta";
            }
            
            //non-singleton group
            else
            {
                fastafile = nonSingletonDir + Utils::getPathSeparator() + Utils::intToString(groupID) + ".fasta";
            }
            
            ofstream fout(fastafile.c_str());
            
            for (vector<string>::const_iterator seq_it = (seqIDs).begin(); seq_it != (seqIDs).end(); seq_it++)
            {
                //cout << ">" << *seq_it << endl;
                fout << ">" << *seq_it << endl; /* write the seqID */
                
                fout << seqMap.at(*seq_it) << endl; /* write the actual sequence */
            }
            
            fout.close();
        }
    }
    
    // align sequences to create msa
    // input: fastafile; output: aligned_fasta
    void align_fasta(const string & align_cmd_prefix, // first part of alignment command
                     string & fastafile, // input: fasta file to align
                     const string & msaDir, // output directory
                     bool clean)
    {
        if (!Utils::isFileEmpty(fastafile))
        {
            // file base with the directory
            string file_base = Utils::getFilebase(fastafile);
            // file base name without directory
            string file_name = Utils::getFilename(file_base);
            // msa file name
            string aligned_fasta = msaDir + "/" + file_name + ".fa";
            
            // mafft command
            string align_cmd = align_cmd_prefix + " " + fastafile + " > " + aligned_fasta;
            
            // system call to align
            const int res = system(align_cmd.c_str());
            if (res != 0 )
                cout << "*** Error: Failed command: " <<  align_cmd << endl;
            
            // clean up unaligned fasta file
            if (clean) {
                remove(fastafile.c_str());
            }
        }

    }
    
    // build hmm for every group. if singleton, build directly;
    // otherwise use mafft do alignment first, then build
    void hmmbuild_jj(const string & hmmbuild_cmd_firstpart, // first part of hmmbuild command
                                                            // hmmbuild --amino --cpu 1 [hmmfile msafile]
                     const string & fastafile, // fasta file with sequences for the group
                     const string & hmmDir, // output directory for all the hmms
                     bool clean)
    {
        // if output directory does not exist, create it
        if (!Utils::isDirectory(hmmDir))
            mkdir(hmmDir.c_str(),0777);
        
        // file base with the directory
        string file_base = Utils::getFilebase(fastafile);
        // file base name without directory
        string file_name = Utils::getFilename(file_base);
        
        
        string hmmfile = hmmDir + "/" + file_name + ".hmm";
        
        //hmmbuild
        string hmmbuild_cmd = hmmbuild_cmd_firstpart + " " + hmmfile + " " + fastafile;
        cout << "hmmbuild command is " << hmmbuild_cmd << endl;
        const int res = system(hmmbuild_cmd.c_str());
        if (res != 0 )
            cout << "*** Error: Failed command: " <<  hmmbuild_cmd << endl;
        
        //remove aligned fasta file
        if (clean) {
            remove(fastafile.c_str());
        }

    }
    
    // split big fasta file into small ones, depending on the number of nodes
    void split_fasta(const map<string, string> & seqMap, const map<int, vector<string> > & lenMap,
                     const map<int, vector<string> > & groupMap, const vector<int> & groupNames,
                     const int & num_nodes,  const string & outputDir)
    {
        
        /********************** generate data for each node, and find chunk information *********************/
        
        /* split big fasta into small fastas */
        vector<vector<int> > chunk_info;
        
        int num_groups, gp_num;
        
        cout << Utils::currentDateTime() << endl;
        
        /* total number of groups */
        num_groups = (int) groupMap.size();
        
        /* number of chunks, cannot be bigger than the total number of groups */
        int num_big_chunks = min(num_nodes, num_groups);
        
        /* number of groups every big chunk has */
        if (num_groups % num_big_chunks == 0) {
            gp_num = num_groups/num_big_chunks;
        }
        else{
            gp_num = num_groups/(num_big_chunks - 1);
            if (num_groups % (num_big_chunks - 1) == 0) {
                num_big_chunks = num_big_chunks - 1;
            }
        }

        
        cout << "   " << num_groups << " total groups: "  << *(groupNames.begin()) << " to " << *(groupNames.rbegin()) << endl;
        cout << "   " << gp_num << " groups every big chunk" << endl;
        cout << "   " << num_nodes << " big chunks" << endl << endl;
        
        Utils::mkdirIfNonExist(outputDir);
        
        // TODO: COULD USE THREADING TO PARALLELIZE THIS PART //
        #pragma omp parallel shared(num_groups, gp_num, num_big_chunks, seqMap, lenMap, groupMap, groupNames, outputDir, chunk_info)
        {
//            int num_threads = omp_get_num_threads();
//            int id = omp_get_thread_num();
//            if(id==0) cout << "number of threads is " << num_threads << endl;
            #pragma omp for
            for (int big_chunk = 0; big_chunk < num_big_chunks; big_chunk++) {
                
                int i, gp_start, gp_end, max_len, min_len, up_len, low_len, gp,chunks,chunk_group, len;
                vector<string> seqs;
                /* starting group for this group, start with 0, and,
                 * ending group for this group, end with num_groups - 1 (included) */
                gp_start = big_chunk * gp_num;
                gp_end = min(num_groups-1, (big_chunk+1) * gp_num - 1);
                //cout << "big chunk " << big_chunk << " start_group: " << gp_start << "; end_group: " << gp_end << endl;
                
                vector<string>::const_iterator it;
                
                /* find the maximum length and upper bound for this node */
                
                min_len = 5000;
                max_len = 0;
                for (i = gp_start; i <= gp_end; i++) {
                    seqs = groupMap.at(groupNames.at(i));
                    for (it = seqs.begin(); it != seqs.end(); it++)
                    {
                        len = (int) seqMap.at(*it).length();
                        min_len = (len < min_len) ? len : min_len;
                        max_len = (len > max_len) ? len : max_len;
                    }

                }
                up_len = (ceil((double)max_len/0.8)> 5000) ? 5000 : ceil((double)max_len/0.8);
                low_len = (floor((double)min_len*0.8) > 30) ? floor((double)min_len*0.8) : 30;
                cout << "   " << big_chunk << " -> min_len: " << min_len << " max_len: " << max_len \
                << " number of groups: " << gp_end - gp_start + 1 << endl;
                
                
                //cout << " Now write the sequences and their groups in the sub file. \n";
                
                map<int, vector<string> >::const_iterator itlow, itup, len_it;
                itlow = lenMap.lower_bound(low_len);
                itup = lenMap.upper_bound(up_len);
                
                //cout << "length that is less than or equal to the min_len: " << itlow->first << endl;
                //cout << "length that is bigger than or equal to the min_len: " << (--itup)->first << endl;
                
                /*************************** write subfiles *****************************/
                string gp_file_base = outputDir + Utils::getPathSeparator() + Utils::intToString(big_chunk);
                string fastafilename = gp_file_base + ".fasta";
                string groupfilename = gp_file_base + ".group";
                
                //cout << fastafilename << " " << groupfilename << endl;
                
                ofstream myfasta(fastafilename.c_str());
                for (len_it = itlow; len_it != itup; len_it++) {
                    for( it = (len_it->second).begin(); it != (len_it->second).end(); it++)
                        myfasta << ">" << *it << endl << seqMap.at(*it) << endl;
                }
                myfasta.close();
                
                ofstream mygroup(groupfilename.c_str());
                for (gp = gp_start; gp <= gp_end; gp++) {
                    for (it = (groupMap.at(groupNames.at(gp))).begin(); it !=  (groupMap.at(groupNames.at(gp))).end(); it++) {
                        mygroup << *it << " " << groupNames.at(gp) << endl;
                    }
                }
                mygroup.close();
                /************************ finish write subfiles **************************/
                //cout << " Finished writing the sequences and their groups in the sub file. \n";
                //cout << endl << endl;
                
                // if the shortest sequence is longer than 60, divide into finer chunks
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

                    cout << "   " << big_chunk << " -> chunks = " << chunks << endl;
                    for (i = 0; i < chunks; i++) {
                        vector<int> this_chunk_info;
                        this_chunk_info.push_back(big_chunk);
                        this_chunk_info.push_back(i*chunk_group); /* start index of this small chunk */
                        this_chunk_info.push_back( min(gps-1, (i+1)*chunk_group - 1) ); /* end index of this small chunk */
                        //cout << this_chunk_info.at(0) << " " << this_chunk_info.at(1) << " " << this_chunk_info.at(2) << endl;
                        #pragma omp critical
                            chunk_info.push_back(this_chunk_info);
                    }
                }
                // if shortest sequence is shorter than 60
                else{
                    vector<int> this_chunk_info;
                    this_chunk_info.push_back(big_chunk);
                    this_chunk_info.push_back(0);
                    this_chunk_info.push_back(gp_end-gp_start);
                    #pragma omp critical
                        chunk_info.push_back(this_chunk_info);
                }
                
            }
        }
        // check
        cout << "   " << chunk_info.size() << " small chunks" << endl;
        string chunk_file = outputDir + Utils::getPathSeparator() + "chunk.info";
        ofstream mychunk(chunk_file.c_str());
        for (vector<vector<int> >::iterator it = chunk_info.begin(); it != chunk_info.end(); it++) {
            mychunk << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << endl;
        }
        mychunk.close();
        
    }



    /********************************************************************/
    /* Parse tab delimited file compatible with the NCBI BLAST m8 and outfmt 6 formats.
     
    Sample record:

    sp|A0A183|LCE6A_HUMAN   tr|G3RZ01|G3RZ01_GORGO  97.5    80      2       0       1       80      1       80      *       *
    sp|A0A183|LCE6A_HUMAN   tr|H2RED1|H2RED1_PANTR  96.2    80      3       0       1       80      1       80      *       *
    sp|A0A183|LCE6A_HUMAN   tr|G1RHA7|G1RHA7_NOMLE  95.0    80      4       0       1       80      1       80      *       *
    sp|A0A183|LCE6A_HUMAN   tr|G7NUF0|G7NUF0_MACFA  93.8    80      5       0       1       80      1       80      *       *
    sp|A0A183|LCE6A_HUMAN   tr|G7MDP4|G7MDP4_MACMU  93.8    80      5       0       1       80      1       80      *       *
    sp|A0A183|LCE6A_HUMAN   tr|F7H175|F7H175_CALJA  86.2    80      11      0       1       80      1       80      *       *
    sp|A0A183|LCE6A_HUMAN   tr|F7HW90|F7HW90_CALJA  85.0    80      12      0       1       80      1       80      *       *

     The fields are: query label, target label, percent identity, alignment length, mismatches, number of gap opens,
     1-based position of start in query, 1-based position of end in query, 1-based position of start in target, position of end in target,
     E-value calculated using Karlin-Altschul statistics, Bit score calculated using Karlin-Altschul statistics
     
     Focus on the third column, find all pairwise identities, and get the mean and median of these.

    */
    void parse_blast6out(const string & filename, // blast6out file name
                      vector<float> & identities, // percentage of identities
                      const bool & clean) //resulting map, remove the blast6out file or not?
    {
        // if file is empty, print message, but don't quit
        if (Utils::isFileEmpty(filename))
        {
            cout << "blast6out file "<< filename << " is empty \n";
        }
        
        //if not either, parse the file
        else
        {
            //myfile -- infile stream
            ifstream myfile(filename.c_str());
            
            char first_char; // first character of the line
            
            while(myfile.get(first_char))
            {
                // ignore lines that start with "#"
                if (first_char == '#') {
                    myfile.ignore(MAX_WIDTH,'\n');
                }
                
                // target sequences
                else
                {
                    myfile.unget(); //put the first character back in the stream
                    
                    string querySeqID, targetSeqID, others;
                    float identity;
                    
                    myfile >> querySeqID >> targetSeqID >> identity; //get identity
                    // ignore the rest of the line after the ID, length, and total score are read
                    myfile.ignore(MAX_WIDTH,'\n');

                    identities.push_back(identity);
                    
                }
            }// end reading file
            
            myfile.close();
            if (clean) {
                remove(filename.c_str());
            }
        }
    }
    
}
