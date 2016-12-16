#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include<math.h> // use sqrt()
#include <vector> // use vectors
#include <algorithm> // use abs()

#include <ctime>

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

// Error detection:
// valgrind ./a.out L030sgGAUSR_e0c100sgGAU_J1_1c0_j2_0c0_G0c80__DBC_J1______map.sd

using namespace std;

int main( int argc, const char* argv[] )
{

	//system("./rm_header.sh");

	std::ofstream OutputFile ;
	OutputFile.open("cat_command.sh", ios::out) ;

	OutputFile << "# /bin/bash " << endl ;
	OutputFile << endl ;
	OutputFile << "# Remove header" << endl ;
	OutputFile << "for seed in $(seq 1 1 445)" << endl ;
	OutputFile << "do" << endl ;
	OutputFile << "sed -i '1d' Tau_Sm1_Sm_$seed.dat" << endl ;
	OutputFile << "done" << endl ;
	OutputFile << endl ;
	OutputFile << "# Join all the files" << endl ;
	OutputFile << "cat " ;
	for (unsigned int k=1; k<445; k++) OutputFile << "Tau_Sm1_Sm_" << k << ".dat " ;
	OutputFile << "> all.dat" << endl ;
	OutputFile << endl ;
	OutputFile << "# Order the file" << endl ;
	OutputFile << "sort -n all.dat > all_f.dat" << endl ;
	OutputFile << endl ;
	OutputFile << "# Render the graph" << endl ;
	OutputFile << "gle -d pdf -cairo smax.gle" << endl ;

	system("chmod +x cat_command.sh");
	system("./cat_command.sh");
	return (0) ;
}
