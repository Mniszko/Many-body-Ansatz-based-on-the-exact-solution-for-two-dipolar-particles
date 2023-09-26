#ifndef MAIN_FUNCTIONS
#define MAIN_FUNCTIONS

#include <array>
#include <vector>

double main_loop(int argc, char* argv[], double g, double g_dip, double length, int nmax, int umax, int mmax, bool WriteWaveFunction = false, double cutoff=0.);

double calculate_full_dispersion(std::vector<std::array<int, 6>>, double*,  bool centerOfMass=false, double cutoff=0.);

//loops over parameters, that call other functions and result in writing output files
//currently changing number of points generated is achieved by changing declarations of these, future plans include further customizing arguments.
int dispersion_loop(int argc, char* argv[], int nmax, int umax, int mmax, bool centerOfMass=false, double cutoff=0.);

int delta_loop(int argc, char* argv[], int nmax, int umax, int mmax, double cutoff=0.);

int length_loop(int argc, char* argv[], int nmax, int umax, int mmax, double cutoff=0.);

int dipole_loop(int argc, char* argv[], int nmax, int umax, int mmax, double cutoff=0.);

#endif