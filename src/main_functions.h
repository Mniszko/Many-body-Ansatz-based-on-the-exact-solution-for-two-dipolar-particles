#ifndef MAIN_FUNCTIONS
#define MAIN_FUNCTIONS

#include <array>
#include <vector>

double main_loop(int argc, char* argv[], double g, double g_dip, double length, int nmax, int umax, int mmax, bool WriteWaveFunction = false);

double calculate_full_dispersion(std::vector<std::array<int, 6>>, double*,  bool centerOfMass=false);

//loops over parameters, that call other functions and result in writing output files
//currently changing number of points generated is achieved by changing declarations of these, future plans include further customizing arguments.
int dispersion_loop(int argc, char* argv[], int nmax, int umax, int mmax, bool centerOfMass=false);

int delta_loop(int argc, char* argv[], int nmax, int umax, int mmax);

int length_loop(int argc, char* argv[], int nmax, int umax, int mmax);

int dipole_loop(int argc, char* argv[], int nmax, int umax, int mmax);

#endif