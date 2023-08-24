#include "./src/integral.h" // stores integral class with all functionalities
#include "./src/variables.h"
#include "./src/flag_parser.h" // used to parse flags
#include "./src/Eigen/Dense"
#include "./src/Eigen/Eigenvalues"
#include <iostream>
#include <cmath>
#include <chrono>
#include <sstream>
#include <exception>
#include <string>
#include <gsl/gsl_errno.h>
#include <vector>
#include <array>
#include <algorithm>

const double pi = 3.14159265358979323846;

double main_loop(int argc, char* argv[], double g, double g_dip, double length, int nmax, int umax, int mmax){
    double epsabs = 1e-8;
    double epsrel = 1e-6;
    size_t limit = 100;
    int nmin = 0;

    int umin = 0;
    
    int mmin = 0;

    int n1 = 0, m1 = 0, u1 = 0, n2 = 0, m2 = 0, u2 = 0, n3 = 0, m3 = 0, u3 = 0, n4 = 0, m4 = 0, u4 = 0;
    
    FlagParser flag_parser(argc, argv);
    flag_parser.parse_flags();

    if (flag_parser.change_lFlag){
        length = flag_parser.change_lValue;
    }
    if (flag_parser.change_gFlag){
        g = flag_parser.change_gValue;
    }    
    if (flag_parser.minmaxFlag){
            nmin = flag_parser.nmin;
            nmax = flag_parser.nmax;
            mmin = flag_parser.mmin;
            mmax = flag_parser.mmax;
            umin = flag_parser.umin;
            umax = flag_parser.umax;
    }


    CompleteIntegral integral(limit,epsabs,epsrel);
    integral.change_length(length);
    std::vector<std::array<int,6>> checker_array;
    integral.change_key(6);

    std::vector<std::array<int, 6>> vector_of_vectors = *generateCombinations(nmin, umin, mmin, nmax, umax, mmax);
    unsigned int number_of_functions = vector_of_vectors.size();
    std::vector<std::array<double,2>> norms;
    for (int i = 0; i<number_of_functions ; ++i){
        std::array<int,6> l1l2 = vector_of_vectors.at(i);
        int n1 = l1l2[0];
        int u1 = l1l2[1];
        int m1 = l1l2[2];
        int n2 = l1l2[3];
        int u2 = l1l2[4];
        int m2 = l1l2[5];
        //std::cout << n1<<' '<<u1<<' '<<m1<<" | "<<n2<<' '<<u1<<' '<<m1 << std::endl;
        //calculating norms in bulk to make performance better (plans for future include finding exact equation that will make the code a bit faster)
        double norm1 = integral.integrate_norm_external(n1, u1, m1);
        double norm2 = integral.integrate_norm_external(n2, u2, m2);
        norms.push_back({1/std::sqrt(norm1),1/std::sqrt(norm2)}); 
    }



    Eigen::MatrixXd main_matrix = Eigen::MatrixXd::Zero(number_of_functions,number_of_functions);

    //first two ifs only consider diagonal elements - updating them for full matrix generation requires some time. Besides with current structure they do nothing.
    if (flag_parser.debugFlag && flag_parser.debugValues.size()==12){ //debugging purpose
        n1 = flag_parser.debugValues.at(0);
        m1 = flag_parser.debugValues.at(1);
        u1 = flag_parser.debugValues.at(2);
        n2 = flag_parser.debugValues.at(3);
        m2 = flag_parser.debugValues.at(4);
        u2 = flag_parser.debugValues.at(5);
        n3 = flag_parser.debugValues.at(0);
        m3 = flag_parser.debugValues.at(1);
        u3 = flag_parser.debugValues.at(2);
        n4 = flag_parser.debugValues.at(3);
        m4 = flag_parser.debugValues.at(4);
        u4 = flag_parser.debugValues.at(5);
        //-------------------------------------------------
        //place debugging code for single tensor element inside


        //-------------------------------------------------

    } else {
        for(int i=0; i<number_of_functions; ++i){
            for(int j=0; j<=i; ++j){

                double result = 0;

                std::array<int,6> l1l2 = vector_of_vectors.at(i);
                std::array<int,6> l3l4 = vector_of_vectors.at(j);

                double norm1 = norms.at(i)[0];
                double norm2 = norms.at(i)[1];
                double norm3 = norms.at(j)[0];
                double norm4 = norms.at(j)[1];

                integral.changeMultiIndex(l1l2[0], l1l2[1], l1l2[2], 1);
                integral.changeMultiIndex(l1l2[3], l1l2[4], l1l2[5], 2);
                integral.changeMultiIndex(l3l4[0], l3l4[1], l3l4[2], 3);
                integral.changeMultiIndex(l3l4[3], l3l4[4], l3l4[5], 4);

                //adding contributions from specific parts of hamiltonian to the result
                
                result += integral.fast_add_over_harmonic(i, j); //computing those energies isn't useful since we know them perfectly well. besides I already wrote the equivalent code using integrals and it slightly worsens performance
                
                //delta potential integration
                result += (integral.integrate_over_delta(g)
                    *integral.get_rev_length()
                    /pi
                    *norm1*norm2*norm3*norm4
                    );

                //dipole potential integration
                /*
                result += (integral.integrate_over_dipole(g_dip)//I'm not sure about that one yet
                    *integral.get_rev_length()
                    *norm1*norm2*norm3*norm4
                    );
                */

                main_matrix(i,j) = result;
                main_matrix(j,i) = result;
            }
        }
        //performing operations on matrix

        //std::cout << "Eigen Matrix:" << std::endl;
        //std::cout << main_matrix << std::endl;

        Eigen::EigenSolver<Eigen::MatrixXd> solver(main_matrix);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Eigenvalue decomposition failed." << std::endl;
            return -1;
        }
        
        Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

        //Finding the lowest eigenvalue
        double minEigenvalue = eigenvalues.minCoeff();
        return minEigenvalue;
    }
    return 0;
}

int length_loop(int argc, char* argv[], int nmax, int umax, int mmax){
    double g = 1.;
    double g_dip = 1.;
    double lengthmax = 10;
    const unsigned int numofsteps = 100;
    std::vector<double> yvalues;
    std::vector<double> xvalues;
    for (double length=1; length < lengthmax ; length+=lengthmax/numofsteps){
        auto start = std::chrono::high_resolution_clock::now();
        double value = main_loop(argc, argv, g, g_dip, length, nmax, umax, mmax);
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_seconds = std::chrono::duration<double>(end - start).count();
        std::cout<< "time elapsed per iteration: "<< elapsed_seconds << " s "<< "iteration number: " << length*numofsteps/lengthmax << '/'<<numofsteps<< std::endl;
        yvalues.push_back(value);
        xvalues.push_back(length);
    }
    write_to_csv("length_output.csv",xvalues,yvalues, nmax, umax, mmax);
    return 0;

    return 0;
}

int delta_loop(int argc, char* argv[], int nmax, int umax, int mmax){
    double g_dip = 1.;
    double length = 1.;
    double gmax = 40.;
    const unsigned int numofsteps = 100;
    std::vector<double> yvalues;
    std::vector<double> xvalues;
    for (double g=0; g < gmax ; g+=gmax/numofsteps){
        auto start = std::chrono::high_resolution_clock::now();
        double value = main_loop(argc, argv, g, g_dip, length, nmax, umax, mmax);
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_seconds = std::chrono::duration<double>(end - start).count();
        std::cout<< "time elapsed per iteration: "<< elapsed_seconds << " s "<< "iteration number: " << g*numofsteps/gmax << '/'<<numofsteps<< std::endl;
        yvalues.push_back(value);
        xvalues.push_back(g);
    }
    write_to_csv("g_output.csv",xvalues,yvalues, nmax, umax, mmax);
    return 0;
}

int dipole_loop(int argc, char* argv[], int nmax, int umax, int mmax){
    std::cout << "NOT YET IMPLEMENTED" <<std::endl;
    return 0;
}

int main(int argc, char* argv[]){
    int nmax = 2;
    int umax = 2;
    int mmax = 1;
    delta_loop(argc, argv, nmax, umax, mmax);
    length_loop(argc, argv, nmax, umax, mmax);
    return 0;
}
