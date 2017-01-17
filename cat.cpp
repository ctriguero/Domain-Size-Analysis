#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h> // use sqrt()
#include <vector> // use vectors
#include <algorithm> // use abs()

#include <ctime>

// Count directory
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>

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
	system("unzip Tau_Smax_data.zip") ;

	// BEGIN Detecting the number of files to process
	system("ls Tau* | wc -l > NUMBER") ;
	ifstream infi ;
	infi.open("NUMBER") ;
	string sLine;
	getline(infi, sLine) ;
    	infi.close() ;
	system("rm -fv NUMBER") ;
	// string to integer:
	int Nfiles = atoi(sLine.c_str()) ;
	cout << "Number of files to be processed: " << Nfiles << endl ;
	int NfilesP1 = Nfiles + 1 ;
	//return (0) ;
	// END Detecting the number of files to process


	std::ofstream OutputFile ;
	OutputFile.open("cat_command.sh", ios::out) ;
	OutputFile << "# /bin/bash " << endl ;
	OutputFile << endl ;
	OutputFile << "# Remove header" << endl ;
	OutputFile << "for seed in $(seq 1 1 "<< Nfiles <<")" << endl ;
	OutputFile << "do" << endl ;
	OutputFile << "sed -i '1d' Tau_Sm1_Sm_$seed.dat" << endl ;
	OutputFile << "done" << endl ;
	OutputFile << endl ;
	OutputFile << "# Join all the files" << endl ;
	OutputFile << "cat " ;
	for (unsigned int k=1; k<NfilesP1; k++) OutputFile << "Tau_Sm1_Sm_" << k << ".dat " ;
	OutputFile << "> all.dat" << endl ;
	OutputFile << endl ;
	OutputFile << "# Order the file" << endl ;
	OutputFile << "sort -n all.dat > all_f.dat" << endl ;
	OutputFile << endl ;
	OutputFile << "# Render the graph" << endl ;
	OutputFile << "gle -d pdf -cairo smax.gle" << endl ;

	system("chmod +x cat_command.sh");
	cout << BOLDRED << "Bash script built & ready to be launched" << RESET << endl ;
	// to execute we use sh ...
	system("sh cat_command.sh");
	system("rm cat_command.sh");
	system("rm all.dat");
	system("rm Tau_*.dat");
	// GLE Graph
	std::ofstream GleFile ;
	GleFile.open("kk.gle", ios::out) ;

	GleFile << "size 13 8" << endl ;
	GleFile << "set texlabels 1" << endl ;
	GleFile << "begin graph" << endl ;
	GleFile << "   scale auto" << endl ;
	GleFile << "   xtitle \"Driving parameter, $\\tau$\" hei 0.5" << endl ;
	GleFile << "   ytitle \"Untransformed space size, $S_{\\rm max}$\" hei 0.5" << endl ;
	GleFile << "   xaxis min -0.5 max -0.37" << endl ;
	GleFile << "   yaxis log min 1 max 11000" << endl ;
	GleFile << "	data \"all_f.dat\" d1 =c1,c22" << endl ;
	GleFile << "	data \"all_f.dat\" d2 =c1,c2" << endl ;
	GleFile << "	data \"all_f.dat\" d3 =c1,c3" << endl ;
	GleFile << "	data \"all_f.dat\" d4 =c1,c4" << endl ;
	GleFile << "	data \"all_f.dat\" d5 =c1,c5" << endl ;
	GleFile << "	data \"all_f.dat\" d6 =c1,c6" << endl ;
	GleFile << "	data \"all_f.dat\" d7 =c1,c7" << endl ;
	GleFile << "	data \"all_f.dat\" d8 =c1,c8" << endl ;
	GleFile << "	data \"all_f.dat\" d9 =c1,c9" << endl ;
	GleFile << "	data \"all_f.dat\" d10 =c1,c10" << endl ;
	GleFile << "	data \"all_f.dat\" d11 =c1,c11" << endl ;
	GleFile << "	d1 deresolve " << Nfiles << " average marker fcircle color black msize 0.05" << endl ; 
	GleFile << "	d2 deresolve " << Nfiles << " average marker fcircle color blue msize 0.05" << endl ; 
	GleFile << "	d3 deresolve " << Nfiles << " average marker fcircle color green msize 0.05" << endl ;    
	GleFile << "	d4 deresolve " << Nfiles << " average marker fcircle color yellow msize 0.05" << endl ;
	GleFile << "	d5 deresolve " << Nfiles << " average marker fcircle color orange msize 0.05" << endl ;    
	GleFile << "	d6 deresolve " << Nfiles << " average marker fcircle color red msize 0.05" << endl ;
	GleFile << "	d7 deresolve " << Nfiles << " average marker fcircle color maroon msize 0.05" << endl ;
	GleFile << "	d8 deresolve " << Nfiles << " average marker fcircle color purple msize 0.05" << endl ;
	GleFile << "	d9 deresolve " << Nfiles << " average marker fcircle color gray msize 0.05" << endl ;             
	GleFile << "	d10 deresolve " << Nfiles << " average marker fcircle color yellowgreen msize 0.05" << endl ;
	GleFile << "	d11 deresolve " << Nfiles << " average marker fcircle color rosybrown msize 0.05" << endl ;
	GleFile << "end graph" << endl ;
	GleFile << "set hei 0.3" << endl ;
	GleFile << "begin key" << endl ;
	GleFile << "nobox" << endl ;
	GleFile << "pos tl" << endl ;
	GleFile << "line color black text \"$S_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "line color blue text \"$S^1_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "line color green text \"$S^2_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "line color yellow text \"$S^3_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "line color orange text \"$S^4_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "line color red text \"$S^5_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "line color maroon text \"$S^6_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "line color purple text \"$S^7_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "line color gray text \"$S^8_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "line color yellowgreen text \"$S^9_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "line color rosybrown text \"$S^{10}_{\\rm max}$\" lwidth 0.1" << endl ;
	GleFile << "end key" << endl ;
	GleFile << "set hei 0.3" << endl ;
	GleFile << "begin key" << endl ;
	GleFile << "nobox" << endl ;
	GleFile << "pos tc" << endl ;
	GleFile << "text \"Cooling, $" << Nfiles << "$ Realizations\"" << endl ;
	GleFile << "text \"$N=100\\times 99$\"" << endl ;
	GleFile << "text \"$\\Delta=0.80$\"" << endl ;
	GleFile << "text \"Average " << Nfiles << "\"" << endl ;
	GleFile << "end key" << endl ;
	GleFile.close() ;

	system("gle -d pdf -cairo kk.gle");
	system("rm -fv kk.gle");
	system("rm -fv all_f.dat");

	return (0) ;
}
