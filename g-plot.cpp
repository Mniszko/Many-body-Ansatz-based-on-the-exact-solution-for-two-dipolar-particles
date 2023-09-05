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
    std::cout << "FINISHED DELTA OF G LOOP" << std::endl;

    //length_loop(argc, argv, nmax, umax, mmax);
    std::cout << "FINISHED DELTA OF LENGTH LOOP" << std::endl;
    
    //dipole_loop(argc, argv, nmax, umax, mmax);
    
    // center of mass dispersion isn't as fast as it could be - for now it computes double integral but could instead be broken into 8 first order integrals which would make the code much faster. If it will be useful I'll implement that
    dispersion_loop(argc, argv, nmax, umax, mmax, true);
    std::cout << "FINISHED CENTER OF MASS LOOP" << std::endl;
    
    //dispersion_loop(argc, argv, nmax, umax, mmax, false);
    std::cout << "FINISHED ONE PARTICLE DISPERSION LOOP" << std::endl;
    return 0;
}

//often main_loop is called to generate already existing output - first it should be checked if the output already exists to make the code around 2 times faster