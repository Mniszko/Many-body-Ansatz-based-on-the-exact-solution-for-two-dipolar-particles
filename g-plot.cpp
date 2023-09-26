#include "./src/integral.h" // stores integral class with all functionalities
#include "./src/variables.h" // additional static functionalities
#include "./src/flag_parser.h" // used to parse flags
#include "./src/main_functions.h" // calculates matrices
#include "./src/Eigen/Dense" // library for matrix manipulation
#include "./src/Eigen/Eigenvalues" // library for matrix manipulation
#include <iostream>




int main(int argc, char* argv[]){

    int nmax = 2;
    int umax = 3;
    int mmax = 3;
    double cutoff = 0.;

    if (argc == 4) { 
        std::cout << "parsing arguments" << std::endl;
        try {
            nmax = std::atoi(argv[1]);
            umax = std::atoi(argv[2]);
            mmax = std::atoi(argv[3]);
            for (int i = 1; i < argc - 3; ++i) {
                argv[i] = argv[i + 3];
            }
            argc -= 3;
            std::cout << "parsed quantum numbers as arguments" << std::endl;
        } catch (...) {
            std::cout << "arguments weren't numbers" << std::endl;
        }

    } else if (argc==2) {
        nmax = -1;
        umax = -1;
        mmax = -1;
        try {
            cutoff = std::atof(argv[1]);
            std::cout << "parsed energy cutoff as argument" << std::endl;
            argc -= 1;
        } catch (...){
            std::cout << "energy wasn't a number" << std::endl;
        }
    }else{
        std::cout << "Warning - no arguments parsed" << std::endl;
    }

    if ((abs(mmax)>umax || umax<0) && nmax!=-1){
        std::cerr << "Parsed impossible quantum numbers " << std::endl;
    }

    std::cout << "maximum quantum numbers: " << nmax << ' ' << umax << ' ' << mmax << std::endl;

    delta_loop(argc, argv, nmax, umax, mmax, cutoff);
    std::cout << "FINISHED DELTA OF G LOOP" << std::endl;

    length_loop(argc, argv, nmax, umax, mmax, cutoff);
    std::cout << "FINISHED DELTA OF LENGTH LOOP" << std::endl;
    
    //dipole_loop(argc, argv, nmax, umax, mmax, cutoff);
    
    // center of mass dispersion isn't as fast as it could be - for now it computes double integral but could instead be broken into 8 first order integrals which would make the code much faster. If it will be useful I'll implement that
    dispersion_loop(argc, argv, nmax, umax, mmax, true, cutoff);
    std::cout << "FINISHED CENTER OF MASS LOOP" << std::endl;
    
    dispersion_loop(argc, argv, nmax, umax, mmax, false, cutoff);
    std::cout << "FINISHED ONE PARTICLE DISPERSION LOOP" << std::endl;
    
    std::cout << "\n\n" << std::endl;
    return 0;





}

//often main_loop is called to generate already existing output - first it should be checked if the output already exists to make the code around 2 times faster