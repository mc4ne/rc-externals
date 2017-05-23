#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
using namespace std;

// Declare constants
static const double pi = 3.14159;
static const double deg2rad = pi/180.;
static const double Mp = 0.9384;  // GeV

// Declare variables
double E0, Ep, EpMin, EpMax, EpStep;
double theta, thetaMin, thetaMax, thetaStep;
double Q2, W2;
int    numEpSteps, numThetaSteps;

ofstream gridFile;
string   gridFileName;
string   gridPath ("input/grids/");

// E0        -> Beam Energy (GeV)
// Ep        -> Scattered electron energy (GeV)
// EpMin     -> Minimum scattered electron energy for grid (GeV)
// EpMax     -> Maximum scattered electron energy for grid (GeV)
// EpStep    -> Steps in scattered electron energy for grid (GeV)
// theta     -> Scattering angle (degrees)
// thetaMin  -> Minimum scattering angle (degrees)
// thetaMax  -> Maximum scattering angle (degrees)
// thetaStep -> Steps in theta for grid (degrees)

// Declare functions
double step_Ep(double EpMin, double EpStep, unsigned int j);
double step_theta(double thetaMin, double thetaStep, unsigned int i);
double calc_Q2(double E0, double Ep, double theta); // Momentum transfer
double calc_W2(double E0, double Ep, double theta); // Invariant mass

int main() {

  cout << "\n Enter the name of the input grid file to create (*.inp): ";
  cin  >> gridFileName;
  cout << "\n Enter the beam energy (GeV): ";
  cin  >> E0;
  cout << "\n Enter the minimum scattered electron energy (GeV): ";
  cin  >> EpMin;
  cout << "\n Enter the maximum scattered electron energy (GeV): ";
  cin  >> EpMax;
  cout << "\n Enter the step size of the scattered electron energy (GeV): ";
  cin  >> EpStep;
  if (EpStep <= 0) {
    cerr << "\n You must enter a value > 0 \n\n Exiting now... \n" << endl;
    exit(1);
  }
  cout << "\n Enter the minimum scattering angle (degrees): ";
  cin  >> thetaMin;
  cout << "\n Enter the maximum scattering angle (degrees): ";
  cin  >> thetaMax;
  cout << "\n Enter the integer step size of the scattering angle (degrees): ";
  cin  >> thetaStep;
  if (thetaStep <= 0) {
    cerr << "\n You must enter a value > 0 \n\n Exiting now... \n" << endl;    
    exit(1);
  }
  cout << "\n";

  gridPath.append(gridFileName.c_str());
  gridFile.open(gridPath.c_str());

  if (!gridFile) {
    cerr << "File could not be opened." << endl;
    exit(1);
  }
  // Maintain historic input structure...
  gridFile << "#\n#\n#\n#\n#\n#\n";

  numEpSteps    = round((EpMax - EpMin) / EpStep); 
  numThetaSteps = round((thetaMax - thetaMin) / thetaStep);
 
  theta = 0.0;
  for (unsigned int i = 0; i <= numThetaSteps; i++) {

    theta = step_theta(thetaMin, thetaStep, i);
    Ep    = 0.0;
    Q2    = 0.0;
    W2    = 0.0;
    for (unsigned int j = 0; j <= numEpSteps; j++) {

      Ep = step_Ep(EpMin, EpStep, j);
      Q2 = calc_Q2(E0, Ep, theta);
      W2 = calc_W2(E0, Ep, theta);
      gridFile << fixed << setprecision(3)  
	       << E0 << "\t" << Ep << "\t" << theta << "\t" << W2 << "\t" << Q2 << endl;
    }  // Ep loop
  }  // theta loop

  gridFile.close();
  return 0;
}  // main


double step_Ep(double EpMin, double EpStep, unsigned int j) {
  return EpMin + double(j) * EpStep;
}

double step_theta(double thetaMin, double thetaStep, unsigned int i) {
  return thetaMin + double(i) * thetaStep;
}

// Momentum transfer
double calc_Q2(double E0, double Ep, double theta) {
  return 4.*E0*Ep*pow(sin(0.5*theta*deg2rad), 2.);
}

// Invariant mass
double calc_W2(double E0, double Ep, double theta) {
  return pow(Mp, 2.) + 2.*Mp*(E0 - Ep) - 4.*E0*Ep*pow(sin(0.5*theta*deg2rad), 2.);
}




