#include "./src/integral.h" // stores integral class with all functionalities
#include "./src/variables.h"
#include "./src/flag_parser.h" // used to parse flags
#include "./src/main_functions.h"
#include "./src/Eigen/Dense"
#include "./src/Eigen/Eigenvalues"
#include <iostream>




int main(int argc, char* argv[]){
    /*
    int nmax = 2;
    int umax = 2;
    int mmax = 1;
    */
    int nmax = 2;
    int umax = 2;
    int mmax = 1;
    //delta_loop(argc, argv, nmax, umax, mmax);
    //length_loop(argc, argv, nmax, umax, mmax);
    //dipole_loop(argc, argv, nmax, umax, mmax);
    dispersion_loop(argc, argv, nmax, umax, mmax, true);
    std::cout << "FINISHED CENTER OF MASS LOOP" << std::endl;
    dispersion_loop(argc, argv, nmax, umax, mmax, false);
    std::cout << "FINISHED ONE PARTICLE DISPERSION LOOP" << std::endl;
    return 0;
}
