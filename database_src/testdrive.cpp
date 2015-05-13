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

	string groupfile = argv[1];
	unordered_map<string, UINT16> seq_gp; 
	map<UINT16, vector<UINT16> > group_len_map;

	/********************** test run of function *********************/

	read_group_sg(groupfile, seq_gp, group_len_map);

	/********************** display work end and time record *********************/

	time_t end_time = time(0); /* ending time in epoch */
	cout<< "Done!" << endl
		<< Utils::currentDateTime() << endl
		<< "Total elapsed time = " << double(end_time - start_time) << " [seconds]" << endl
		<< "============================================================================"
		<< endl << endl;


	return 0;
}
