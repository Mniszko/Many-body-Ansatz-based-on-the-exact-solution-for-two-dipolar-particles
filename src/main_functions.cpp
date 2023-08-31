#include "integral.h" // stores integral class with all functionalities
#include "variables.h"
#include "flag_parser.h" // used to parse flags
#include "main_functions.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
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
#include <fstream>


const double pi = 3.14159265358979323846;

double main_loop(int argc, char* argv[], double g, double g_dip, double length, int nmax, int umax, int mmax, bool WriteWaveFunction){
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
        //saving matrix(debugging)
        std::ofstream file("./testy/sample_matrix.csv");
        if (file.is_open()) {
            for (int row = 0; row < main_matrix.rows(); ++row) {
                for (int col = 0; col < main_matrix.cols(); ++col) {
                    file << main_matrix(row, col);
                    if (col + 1 < main_matrix.cols())
                        file << ", ";  // Add comma for all but last element
                }
                file << "\n";  // Newline for next row
            }
        } else {
            std::cerr << "Failed to open the file: " << "./testy/sample_matrix.csv" << std::endl;
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
        double minEigenvalue = eigenvalues.minCoeff();
        
        if (WriteWaveFunction){
            //Finding the lowest eigenvalue
            int minIndex = 0;

            for (int i = 0; i < eigenvalues.size(); ++i) {
                if (eigenvalues[i] == minEigenvalue) {
                    minIndex = i;
                }
            }
            int size = main_matrix.rows();
            double* stdArray = new double[size];


            auto minEigenvector = solver.eigenvectors();
            for (int i = 0; i < size; ++i) {
                stdArray[i] = minEigenvector.col(minIndex).real()(i);
            }
            write_vectors_to_csv(vector_of_vectors, stdArray, number_of_functions, "./vector_outputs/lowest_wavefunction"+std::to_string(length)+".csv");
            delete[] stdArray;
        }

        return minEigenvalue;
    }
    return 0;
}

//calculate_full_dispersion is very simmilar to main_loop but lacks option to debug using flags
double calculate_full_dispersion(std::vector<std::array<int, 6>> vector_of_vectors, double* array, bool centerOfMass){
    int number_of_functions = vector_of_vectors.size();
    double epsabs = 1e-8; //maybe making number if recursions higher would be a good idea
    double epsrel = 1e-6;
    size_t limit = 100;

    CompleteIntegral integral(limit,epsabs,epsrel);

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
    
    Eigen::MatrixXd r_matrix = Eigen::MatrixXd::Zero(number_of_functions,number_of_functions);
    Eigen::MatrixXd r2_matrix = Eigen::MatrixXd::Zero(number_of_functions,number_of_functions);
    for(int i=0; i<number_of_functions; ++i){
            for(int j=0; j<number_of_functions; ++j){

                double result = 0;

                std::array<int,6> l1l2 = vector_of_vectors.at(i);
                std::array<int,6> l3l4 = vector_of_vectors.at(j);

                double c1 = std::abs(array[i]);
                double c2 = std::abs(array[j]);

                double norm1 = norms.at(i)[0];
                double norm2 = norms.at(i)[1];
                double norm3 = norms.at(j)[0];
                double norm4 = norms.at(j)[1];

                integral.changeMultiIndex(l1l2[0], l1l2[1], l1l2[2], 1);
                integral.changeMultiIndex(l1l2[3], l1l2[4], l1l2[5], 2);
                integral.changeMultiIndex(l3l4[0], l3l4[1], l3l4[2], 3);
                integral.changeMultiIndex(l3l4[3], l3l4[4], l3l4[5], 4);
                //if we are integrating over only one-particle coordinate the matrix of components becomes diagonal - otherwise its rather improbable
                if(i==j && !centerOfMass){
                    r_matrix(i,j) = integral.integrate_delta_dispersion_helper_r(c1,c2, centerOfMass)*norm1*norm2*norm3*norm4;//definition of norm inside this function validates this form
                    r_matrix(j,i) = r_matrix(i,j);
                    r2_matrix(i,j) = integral.integrate_delta_dispersion_helper_r2(c1,c2, centerOfMass)*norm1*norm2*norm3*norm4;
                    r2_matrix(j,i) = r2_matrix(i,j);
                } else if (i!=j && !centerOfMass) {
                    r_matrix(i,j) = 0;
                    r_matrix(j,i) = 0;
                    r2_matrix(i,j) = 0;
                    r2_matrix(j,i) = 0;
                } else if (centerOfMass) {
                    r_matrix(i,j) = integral.integrate_delta_dispersion_helper_r(c1,c2, centerOfMass)*norm1*norm2*norm3*norm4;
                    r_matrix(j,i) = r_matrix(i,j);
                    r2_matrix(i,j) = integral.integrate_delta_dispersion_helper_r2(c1,c2, centerOfMass)*norm1*norm2*norm3*norm4;
                    r2_matrix(j,i) = r2_matrix(i,j);
                } else {
                    std::cerr << "Matrix estimation error" << std::endl;
                }
/*
                if (r_matrix(i,j)>=1 || r2_matrix(i,j)>=1){
                    std::cout << "r value: " << r_matrix(i,j) << " r^2 value: " << r2_matrix(i,j) << " c1, c2: " << c1 <<' '<< c2  << std::endl;
                    std::cout << "frst n1,u1,m1,n2,u2,m2:";
                    for (int k=0; k<6 ; ++k){
                        std::cout << l1l2[k]; 
                    }
                    std::cout << std::endl;
                    std::cout << "scnd n1,u1,m1,n2,u2,m2:";
                    for (int k=0; k<6 ; ++k){
                        std::cout << l3l4[k]; 
                    }
                    std::cout << std::endl;
                }
*/
            }
        }
    double rsum = r_matrix.sum();
    double r2sum = r2_matrix.sum();

    //rsum=r_matrix(0,0);
    //r2sum = r2_matrix(0,0);
    std::cout << "r2sum = " <<r2sum << std::endl;
    std::cout << "rsum = " << rsum << std::endl;
    double radial_dispersion = sqrt(r2sum-rsum*rsum);
    return radial_dispersion;
}

int dispersion_loop(int argc, char* argv[], int nmax, int umax, int mmax, bool centerOfMass){
    double g = 1.;
    double g_dip = 1.;
    double lengthmax = 20;
    double length0 = 0.1;
    const unsigned int numofsteps = 100;
    std::vector<double> yvalues;
    std::vector<double> xvalues;

    int nmin = 0;
    int umin = 0;
    int mmin = 0;
    //int n1 = 0, m1 = 0, u1 = 0, n2 = 0, m2 = 0, u2 = 0, n3 = 0, m3 = 0, u3 = 0, n4 = 0, m4 = 0, u4 = 0;
    FlagParser flag_parser(argc, argv);
    flag_parser.parse_flags();
    std::vector<std::array<int, 6>> vector_of_vectors = *generateCombinations(nmin, umin, mmin, nmax, umax, mmax);
    unsigned int number_of_functions = vector_of_vectors.size();
    deleteFilesWithPrefix();
    for (double length=length0; length < lengthmax ; length+=(lengthmax-length0)/numofsteps){
        main_loop(argc, argv, g, g_dip, length, nmax, umax, mmax, true);
        double* array = readLastColumn("./vector_outputs/lowest_wavefunction"+std::to_string(length)+".csv", number_of_functions);
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "-----------------------------------------" << std::endl;
        double value = calculate_full_dispersion(vector_of_vectors, array, centerOfMass);
        std::cout << "result: " << value << "\tfile read:" << length <<std::endl;
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_seconds = std::chrono::duration<double>(end - start).count();
        std::cout<< "time elapsed per iteration: "<< elapsed_seconds << " s "<< "iteration number: " << int(length*numofsteps/lengthmax) << '/'<<numofsteps<< std::endl;
        yvalues.push_back(value);
        xvalues.push_back(length);
    }
    char nameAddon[3] = {'C','O','M'};
    if(!centerOfMass){
        nameAddon[2] = 'G';
    }
    std::string filename = "dispersion_output";
    filename.append(nameAddon);
    filename.append(".csv");
    write_to_csv(filename,xvalues,yvalues, nmax, umax, mmax);
    return 0;
}

int delta_loop(int argc, char* argv[], int nmax, int umax, int mmax){
    double g_dip = 1.;
    double length = 1.;
    double gmax = 40.;
    const unsigned int numofsteps = 20;
    std::vector<double> yvalues;
    std::vector<double> xvalues;
    deleteFilesWithPrefix();
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

int length_loop(int argc, char* argv[], int nmax, int umax, int mmax){
    double g = 1.;
    double g_dip = 10000.;
    double lengthmax = 100;
    const unsigned int numofsteps = 20;
    std::vector<double> yvalues;
    std::vector<double> xvalues;
    deleteFilesWithPrefix();
    for (double length=1; length < lengthmax ; length+=lengthmax/numofsteps){
        auto start = std::chrono::high_resolution_clock::now();
        double value = main_loop(argc, argv, g, g_dip, length, nmax, umax, mmax,true);
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_seconds = std::chrono::duration<double>(end - start).count();
        std::cout<< "time elapsed per iteration: "<< elapsed_seconds << " s "<< "iteration number: " << length*numofsteps/lengthmax << '/'<<numofsteps<< std::endl;
        yvalues.push_back(value);
        xvalues.push_back(length);
    }
    write_to_csv("length_output.csv",xvalues,yvalues, nmax, umax, mmax);
    return 0;
}

int dipole_loop(int argc, char* argv[], int nmax, int umax, int mmax){
    double g = 1.;
    double g_dip = 1.;
    double lengthmax = 10;
    const unsigned int numofsteps = 20;
    std::vector<double> yvalues;
    std::vector<double> xvalues;
    deleteFilesWithPrefix();//cleans vector outputs folder
    for (double length=0.01; length < lengthmax ; length+=lengthmax/numofsteps){
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
}