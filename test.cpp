#include "src/integral.h" // stores integral class with all functionalities
#include "src/variables.h"
#include "src/flag_parser.h" // used to parse flags
#include "src/main_functions.h"
#include "src/Eigen/Dense"
#include "src/Eigen/Eigenvalues"
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
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <fstream>

double h=1.;

double fast_power(double a, int b){
    double result = 1;
    for(int i=0; i<b ; ++i){
        result *= a;
    }
    return result;
}

int main(){
    int nmin=0;
    int umin=0;
    int mmin=0;
    int nmax=1;
    int umax=1;
    int mmax=1;

    std::vector<std::array<double,2>> norms;

    double epsabs = 1e-8;
    double epsrel = 1e-6;
    size_t limit = 100;

    double c1=1,c2=1;

    CompleteIntegral integral(limit,epsabs,epsrel);
    std::vector<std::array<int,6>> combination = *generateCombinations(nmin,umin,mmin,nmax,umax,mmax);
    unsigned int number_of_functions = combination.size();
    Eigen::MatrixXd changer(number_of_functions,number_of_functions);
    integral.change_length(1);
    for (int i=0 ; i < number_of_functions ; ++i){
        for (int j=0 ; j <= i ; ++j){
            std::array<int,6> l1l2 = combination.at(i);
            std::array<int,6> l3l4 = combination.at(j);

            integral.changeMultiIndex(l1l2[0], l1l2[1], l1l2[2], 1);
            integral.changeMultiIndex(l1l2[3], l1l2[4], l1l2[5], 2);
            integral.changeMultiIndex(l3l4[0], l3l4[1], l3l4[2], 3);
            integral.changeMultiIndex(l3l4[3], l3l4[4], l3l4[5], 4);
            double element = integral.fast_add_over_harmonic(i,j);
            element += integral.integrate_over_delta(1);
            changer(i,j) = element;
            changer(j,i) = element;
        }
    }


    std::ofstream file("./testy/sample_matrix.csv");
    if (file.is_open()) {
        for (int row = 0; row < changer.rows(); ++row) {
            for (int col = 0; col < changer.cols(); ++col) {
                file << changer(row, col);
                if (col + 1 < changer.cols())
                    file << ", ";  // Add comma for all but last element
            }
            file << "\n";  // Newline for next row
        }
    } else {
        std::cerr << "Failed to open the file: " << "./testy/sample_matrix.csv" << std::endl;
    }

    for (int i=0 ; i < number_of_functions ; ++i){
        std::cout << combination.at(i)[0] << combination.at(i)[1] << combination.at(i)[2] << combination.at(i)[3] << combination.at(i)[4] << combination.at(i)[5] << std::endl;
    }
    std::cout << changer;
    return 1;
}

/*
int main(){
    std::vector<std::array<double,2>> norms;

    double epsabs = 1e-8;
    double epsrel = 1e-6;
    size_t limit = 100;

    double c1=1,c2=1;

    CompleteIntegral integral(limit,epsabs,epsrel);

    int n1=0;
    int u1=0;
    int m1=0;

    int n2=0;
    int u2=0;
    int m2=0;

    int n3=0;
    int u3=0;
    int m3=0;

    int n4=0;
    int u4=0;
    int m4=0;

    double norm1 = integral.integrate_norm_external(n1, u1, m1);
    double norm2 = integral.integrate_norm_external(n2, u2, m2);
    double norm3 = integral.integrate_norm_external(n3, u3, m3);
    double norm4 = integral.integrate_norm_external(n4, u4, m4);

    integral.changeMultiIndex(n1, u1, m1, 1);
    integral.changeMultiIndex(n2, u2, m2, 2);
    integral.changeMultiIndex(n3, u3, m3, 3);
    integral.changeMultiIndex(n4, u4, m4, 4);


    double r_value = integral.integrate_delta_dispersion_helper_r(c1,c2)/sqrt(norm1)/sqrt(norm2)/sqrt(norm3)/sqrt(norm4);
    double r2_value = integral.integrate_delta_dispersion_helper_r2(c1,c2)/sqrt(norm1)/sqrt(norm2)/sqrt(norm3)/sqrt(norm4);

    double rho1=1.;
    
    std::cout << "r value:\t" << r_value << std::endl;
    std::cout << "r2 value:\t" << r2_value << std::endl;
    std::cout << "mean value:\t" << sqrt(r2_value-r_value*r_value) << std::endl;

    double value_r; 
    return 1;
}
*/
/*
int main(){
    int length_of_vectors = 3;
    std::vector<std::array<int,6>> vector_of_vectors;
    for (int i = 0 ; i < length_of_vectors ; ++i){
        vector_of_vectors.push_back({0,0,0,0,0,i});
    }
    Eigen::MatrixXd main_matrix(length_of_vectors,length_of_vectors);
    main_matrix << 1, 2, 3,
     4, 5, 6,
     7, 8, 9;

    Eigen::EigenSolver<Eigen::MatrixXd> solver(main_matrix);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Eigenvalue decomposition failed." << std::endl;
            return -1;
        }
    Eigen::VectorXd eigenvalues = solver.eigenvalues().real();
    double minEigenvalue = eigenvalues.minCoeff();

    if (true){
        //Finding the lowest eigenvalue
        int minIndex = 0;

        for (int i = 0; i < eigenvalues.size(); ++i) {
            std::cout<< eigenvalues[i] << std::endl;
            if (eigenvalues[i] == minEigenvalue) {
                minIndex = i;
            }
        }
        int size = main_matrix.rows();
        double* stdArray = new double[size];


        auto minEigenvector = solver.eigenvectors();
        std::cout << minEigenvector << std::endl;
        for (int i = 0; i < size; ++i) {
            stdArray[i] = minEigenvector.col(minIndex).real()(i);
        }
        write_vectors_to_csv(vector_of_vectors, stdArray, length_of_vectors, "./testy/test.csv");
        delete[] stdArray;
    }
    return 1;
}
*/