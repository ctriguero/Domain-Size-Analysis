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

	// Parameters for the calculation:
	//------------------------------------------------------------------
	//------------------------------------------------------------------
	int Nc=5 ; // Number of stored cluster sizes
	string SaveS0             = "no" ;    // yes to store positions of s=0 sites
	string Save_S_cooling_ppm = "yes" ;    // yes to directly store S-field images
	string SaveSp1            = "no" ;    // yes to store positions of s=+1 sites
	string SaveSm1            = "no" ;    // yes to store positions of s=-1 sites
	//------------------------------------------------------------------
	//------------------------------------------------------------------
	// Parameters for the calculation:

	int start_s=clock();
	
	cout << endl ;
	cout << endl ;
	

	string line ;
    
	// Variables to be read
	 int xpos, ypos, time, Naval ;
	int new_s, new_d ;
	double PreTau, h, e, e1, e2;

	// BEGIN pipe in file to analyze
	std::ostringstream ss;
	ss << argv[1] << endl ;
	std::ifstream InputFile;
	cout << "\t-> Analyzing the file: " << BOLDRED << argv[1] << RESET << endl ;
//	cout << "\t-> Analyzing the file: " << BOLDRED << ss.str().c_str() << RESET << endl ;
	InputFile.open(argv[1], ios::in) ;
//	InputFile.open("L030sgGAUSR_e0c100sgGAU_J1_1c0_j2_0c0_G0c80__DBC_J1______map.sd", ios::in) ;
	// END pipe in file to analyze

	
	// BEGIN Number of transition states
	 int block_num=0 ;
	while (getline(InputFile, line))
	{
		if (line.length() > 2)
		{		
			istringstream iss1(line); //procesar primera línea (leemos new_s)
			iss1 >> xpos >> ypos >> time >> Naval >> new_s >> new_d ;
			if ( new_s != 0 ) block_num++ ;
			getline(InputFile, line) ; //procesar segunda línea (sin hacer nada)			
			getline(InputFile, line) ; //procesar tercera línea (sin hacer nada)
		}
	   
	}
	InputFile.close() ;
	 int L=(1+sqrt(1+4*block_num))/2 ;
	cout << "\t-> The size of the system is: " << BOLDRED << L << " X " << L-1 << RESET << endl ;
	// END Number of transition states
	

	// BEGIN Initial s=0 sites
	size_t InitialSitesNumber = L*(L-1) ;
	std::vector<int> X(InitialSitesNumber), Y(InitialSitesNumber) ;
	 int counter=0 ;
	for (xpos=0; xpos<L; xpos++)
	{
		for (ypos=0; ypos<L-1; ypos++) 
		{
			X[counter]=xpos ;
			Y[counter]=ypos ;
			counter++ ;
		}
	}
	std::cout << "\t-> Initial s=0 sites detected: " << BOLDRED << X.size() << RESET << endl ;
	// END Initial s=0 sites


	// BEGIN Initial s=+1 sites empty
	std::vector<int> Xp, Yp ;
	// END Initial s=+1 sites empty

	// BEGIN Initial s=-1 sites empty
	std::vector<int> Xm, Ym ;
	// END Initial s=+1 sites empty


	// BEGIN Store cluster information
	ofstream OutputSizes ;
	OutputSizes.open("Tau_Sm1_Sm.dat", ios::out | ios::trunc) ;
	OutputSizes << "!Tau" << "\t\t" ;
	for ( int j=1; j<Nc+1; j++) OutputSizes << "Smax" << j << "\t" ;
	OutputSizes << "Smax" << "\t" << endl ;
	// END Store cluster information

	// BEGIN Read loop
//	InputFile.open("L030sgGAUSR_e0c100sgGAU_J1_1c0_j2_0c0_G0c80__DBC_J1______map.sd", ios::in) ;
	InputFile.open(argv[1], ios::in) ;
	
	int FrameCounter=0 ;
	double PreTauOld=-0.0000001 ;
	block_num=0 ;
	while (getline(InputFile, line))
	{
		if (line.length() > 2)
		{
			istringstream iss1(line); //procesar primera línea
			iss1 >> xpos >> ypos >> time >> Naval >> new_s >> new_d ;
			getline(InputFile, line); //procesar segunda línea
			istringstream iss2(line);
			iss2 >> PreTau >> h >> e ;			
			getline(InputFile, line); //procesar tercera línea (sin hacer nada)


			

			// BEGIN For cooling
			if ( new_s != 0 )
			{
				block_num++;
				cout << BOLDBLUE << "\t-> Block number: " << RESET << block_num << endl ;

				if ( PreTau != PreTauOld )
				{
					// BEGIN Clustering
					std::vector<int> C(X.size()) ; 
					for ( int j=0; j<X.size(); j++) C[j]=j ;
					std::vector<int> auxiliar(X.size()+1,0) ;
					for ( int j=0; j<X.size()-1; j++)
					{
						for ( int k=j+1; k<X.size(); k++)
						{
							// Clustering now using Manhattan metric d=|x2-x1|+|y2-y1|+|z3-z1|
							if ( abs(X[k]-X[j]) + abs(Y[k]-Y[j]) == 1 )
							{		
								unsigned int OldClusterSize=0 ;
								for ( int m=0; m<X.size(); m++)
								{
									if ( C[k] == C[m] )
									{
										OldClusterSize++ ;
										auxiliar[OldClusterSize]=m ;
									}
								}
								for ( int m=1; m<OldClusterSize+1; m++)
								{
									C[ auxiliar[m] ]=C[j] ;
								}
							}
						}
					}
					// END Clustering

					
					// BEGIN Clustering recognition
					for ( int j=0; j<X.size(); j++) auxiliar[j]=0 ;
					 int Number_clusters=0 ;
					for ( int j=0; j<X.size(); j++)
					{
						int skip=0 ;
						for ( int k=0; k<X.size(); k++)
						{
							if ( ( C[k] == j ) && ( skip == 0 ) )
							{
								auxiliar[Number_clusters]=j ;
								Number_clusters++ ;
								skip=1 ;
							}
						}
					}
					cout << "\t-> Detected " << BOLDRED << Number_clusters << " clusters" << RESET << endl ;
					// END Clustering recognition
					
					// BEGIN Store clusters
					std::vector<int> size(Number_clusters,0) ;
					 int Cluster_Index=0 ;
					for ( int j=0; j<Number_clusters; j++)
					{
						 int cluster_atoms=0 ;
						for ( int k=0; k<X.size(); k++)
						{
							if ( C[k] == auxiliar[j] )
							{
								cluster_atoms++ ;
							}
						}
						size[j]=cluster_atoms ;
					}
					std::sort(size.rbegin(), size.rend());
					if ( size.size() == 1)
					{
						cout << "\t   * The cluster contains " << size[0] << " nodes." << endl ;						
					}
					else
					{
						cout << "\t   * Largest cluster contains " << size[0] << " nodes." << endl ;
						cout << "\t   * Smallest cluster contains " << size[Number_clusters-1] << " nodes." << endl ;
					}


					// Here we store the first Nc largest clusters
					double TauOld=PreTauOld*100.0 ;
					OutputSizes << TauOld << "\t" ;
					for ( int j=1; j<Nc; j++)
					{
						if ( size.size() == j)
						{
							for ( int k=0; k<j; k++) OutputSizes << size[k] << "\t" ;
							for ( int k=j; k<Nc; k++) OutputSizes << "0" << "\t" ;
						}
					}
					if ( size.size() >= Nc) 
					{
						for ( int k=0; k<Nc; k++) OutputSizes << size[k] << "\t" ;
					}
					OutputSizes << X.size() << endl ;					
					// END Store clusters

					FrameCounter++ ;

////					// BEGIN Store s=0 maps
					if ( SaveS0 == "yes" )
					{
						std::ostringstream name ;
						
						if ( FrameCounter <= 9) name << "Frame_no_000" << FrameCounter << "_map_s0.xy" ;
						if ( ( FrameCounter >= 10 ) && ( FrameCounter <= 99 ) ) name << "Frame_no_00" << FrameCounter << "_map_s0.xy" ;
						if ( ( FrameCounter >= 100 ) && ( FrameCounter <= 999) ) name << "Frame_no_0" << FrameCounter << "_map_s0.xy" ;
						if ( ( FrameCounter >= 1000 ) && ( FrameCounter <= 9999) ) name << "Frame_no_" << FrameCounter << "_map_s0.xy" ;

						cout << YELLOW << "\t-> Storing map s=0: " << name.str().c_str() << RESET << endl ;
				                std::ofstream OutMapS0 ;
				                OutMapS0.open(name.str().c_str(), ios::out | ios::trunc) ; // .bin  | ios::binary
						OutMapS0 << "#Positions of s=0 sites @ Tau= " << TauOld << endl ;
						OutMapS0 << "#X\tY" << endl ;
						for ( int k=0; k<X.size(); k++) OutMapS0 << X[k] << "\t" << Y[k] << endl ;
						OutMapS0.close() ;
					}
					// END Store s=0 maps


					// BEGIN Store S-field ppm images
					if ( Save_S_cooling_ppm == "yes" )
					{
						std::ostringstream nameppm ;
						
						if ( FrameCounter <= 9) nameppm << "Frame_no_000" << FrameCounter << ".ppm" ;
						if ( ( FrameCounter >= 10 ) && ( FrameCounter <= 99 ) ) nameppm << "Frame_no_00" << FrameCounter << ".ppm" ;
						if ( ( FrameCounter >= 100 ) && ( FrameCounter <= 999) ) nameppm << "Frame_no_0" << FrameCounter << ".ppm" ;
						if ( ( FrameCounter >= 1000 ) && ( FrameCounter <= 9999) ) nameppm << "Frame_no_" << FrameCounter << ".ppm" ;

						cout << YELLOW << "\t-> Storing S-field ppm image: " << nameppm.str().c_str() << RESET << endl ;
				                std::ofstream OutMapSppm ;
				                OutMapSppm.open(nameppm.str().c_str(), ios::out | ios::trunc) ;
						OutMapSppm << "P3" << endl ;
						OutMapSppm << L << " " << L-1 << endl ;
						OutMapSppm << "255" << endl ; // RGB triplets colors: 0-255
						OutMapSppm << "# S-field for cooling" << endl ;
						
						int kkonter=0 ;
						for ( int ii=0; ii<L; ii++)
						{
							for ( int jj=0; jj<L-1; jj++)
							{											
								for ( int kk=0; kk<X.size(); kk++)  if ( ( ii == X[kk] )  && ( jj == Y[kk] ))  OutMapSppm << "0   82  33 " << endl ; // 0 82 33 or 26 102 46
								for ( int kk=0; kk<Xp.size(); kk++) if ( ( ii == Xp[kk] ) && ( jj == Yp[kk] )) OutMapSppm << "153 0   0  " << endl ; 
								for ( int kk=0; kk<Xm.size(); kk++) if ( ( ii == Xm[kk] ) && ( jj == Ym[kk] )) OutMapSppm << "70  130 180" << endl ; 
								kkonter++ ;
							}
						}
						//  ffmpeg -i Frame_no_%04d_s0.ppm -c:v libx264 -r 30 out.mp4
						OutMapSppm.close() ;
						cout << "kkonter= " << kkonter << endl ;
					}
					// END Store S-field images


					// BEGIN Store s=+1 maps
					if ( SaveSp1 == "yes" )
					{
						std::ostringstream name ;

			        	        if ( FrameCounter <= 9) name << "Frame_no_000" << FrameCounter << "_map_sp1.xy" ;
						if ( ( FrameCounter >= 10 ) && ( FrameCounter <= 99 ) ) name << "Frame_no_00" << FrameCounter << "_map_sp1.xy" ;
						if ( ( FrameCounter >= 100 ) && ( FrameCounter <= 999) ) name << "Frame_no_0" << FrameCounter << "_map_sp1.xy" ;
						if ( ( FrameCounter >= 1000 ) && ( FrameCounter <= 9999) ) name << "Frame_no_" << FrameCounter << "_map_sp1.xy" ;

						cout << YELLOW << "\t-> Storing map s=+1: " << name.str().c_str() << RESET << endl ;
				                std::ofstream OutMapSp1 ;
				                OutMapSp1.open(name.str().c_str(), ios::out | ios::trunc) ;
						OutMapSp1 << "#Positions of s=+1 sites @ Tau= " << TauOld << endl ;
						OutMapSp1 << "#X\tY" << endl ;
						for ( int k=0; k<Xp.size(); k++) OutMapSp1 << Xp[k] << "\t" << Yp[k] << endl ;
						OutMapSp1.close() ;
					}
					// END Store s=+1 maps


					// BEGIN Store s=-1 maps
					if ( SaveSm1 == "yes" )
					{
						std::ostringstream name ;
			        	        
						if ( FrameCounter <= 9) name << "Frame_no_000" << FrameCounter << "_map_sm1.xy" ;
						if ( ( FrameCounter >= 10 ) && ( FrameCounter <= 99 ) ) name << "Frame_no_00" << FrameCounter << "_map_sm1.xy" ;
						if ( ( FrameCounter >= 100 ) && ( FrameCounter <= 999) ) name << "Frame_no_0" << FrameCounter << "_map_sm1.xy" ;
						if ( ( FrameCounter >= 1000 ) && ( FrameCounter <= 9999) ) name << "Frame_no_" << FrameCounter << "_map_sm1.xy" ;

						cout << YELLOW << "\t-> Storing map s=-1: " << name.str().c_str() << RESET << endl ;
				                std::ofstream OutMapSm1 ;
				                OutMapSm1.open(name.str().c_str(), ios::out | ios::trunc) ; 
						OutMapSm1 << "#Positions of s=-1 sites @ Tau= " << TauOld << endl ;
						OutMapSm1 << "#X\tY" << endl ;
						for ( int k=0; k<Xm.size(); k++) OutMapSm1 << Xm[k] << "\t" << Ym[k] << endl ;
						OutMapSm1.close() ;
					}
					// END Store s=-1 maps

				}
			}
			PreTauOld=PreTau ;
			// END For cooling
			
			// BEGIN Remove a site
			for ( int i=0; i<X.size(); i++) 
			{
				if ( ( X[i] == xpos ) && ( Y[i] == ypos ) )
				{
					cout << "Site removed is: " << X[i] << " " << Y[i] << endl ;
					X.erase (X.begin()+i) ;
					Y.erase (Y.begin()+i) ;
				}
			}
			// END Remove a site


			// BEGIN add site to vectors s=+1
			if ( new_s == 1 )
			{
				Xp.push_back (xpos) ;
				Yp.push_back (ypos) ;
			}
			// END add site to vectors s=+1


			// BEGIN add site to vectors s=-1
			if ( new_s == -1 )
			{
				Xm.push_back (xpos) ;
				Ym.push_back (ypos) ;
			}
			// END add site to vectors s=-1
		
			
		}
	   
	}
	InputFile.close() ;
	// END Read loop
	cout << "Size should be zero now and is: " << X.size () << endl ;
	
	int stop_s=clock();
	cout << endl ;
	cout << "Total execution time: " << BOLDBLACK << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << RESET << endl;
	
	return (0) ;
}
