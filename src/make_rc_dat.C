#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
using namespace std;

// Declare variables
ifstream externalsOut;
ofstream radCorrFile;
string   externalsOutName, radCorrFileName;
string externalsOutPath ("output/externals/");
string radCorrOutPath   ("output/rad-corr-data/");

double E0, Ep, theta, sigmaBorn, sigmaRad, corrFactor, dummy;

int main() {

  // Define the output path
  // Read in the input file
  cout << "\n Enter the name of the externals output file to read (*.out): ";
  cin  >> externalsOutName;
  externalsOutPath.append(externalsOutName.c_str());
  externalsOut.open(externalsOutPath.c_str());
  if (!externalsOut) {
    cerr << "\n Error reading externals output file! \n" << endl;
    exit(1);
  }
  // Open the output data file
  cout << "\n Enter the name of the data file to create (*.dat): ";
  cin  >> radCorrFileName;
  radCorrOutPath.append(radCorrFileName.c_str());
  radCorrFile.open(radCorrOutPath.c_str());
  if (!radCorrFile) {
    cerr << "\n Output data file could not be opened! \n" << endl;
    exit(1);
  }

  dummy = 0.0; E0 = 0.0; Ep = 0.0; theta = 0.0; 
  sigmaBorn = 0.0; sigmaRad = 0.0; corrFactor = 0.0;
   while (!externalsOut.eof()) {
  // while (externalsOut >> E0 >> Ep >> theta >> dummy >> dummy >> sigmaBorn >> sigmaRad
  // 	 >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
  // 	 >> dummy >> dummy >> dummy) {
    
     externalsOut >> E0 >> Ep >> theta >> dummy >> dummy >> sigmaBorn >> sigmaRad
  		 >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
  		 >> dummy >> dummy >> dummy;
    if (externalsOut.eof()) break;
    corrFactor = sigmaBorn / sigmaRad;
    radCorrFile << fixed << setprecision(3) << E0 << "\t" << Ep << "\t" << theta << "\t" 
  	     << setprecision(5) << sigmaBorn << "\t" << corrFactor << endl;
  } 

  externalsOut.close(); 
  radCorrFile.close();
  return 0;

}  // main
