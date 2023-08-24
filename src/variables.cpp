#include "variables.h"

#include <fstream>
#include <vector>
#include <string>
#include <array>
#include <iostream>
#include <boost/program_options.hpp>
#include <stdexcept>


void functions::write_to_file(std::string filename, const std::vector<std::string>& parameters, const std::vector<std::vector<double>>& values, std::string extrainfo) {
        std::ofstream file(filename); // Open file

        // Write parameters to the file
        file << extrainfo << '\n';
        for (size_t i = 0; i < parameters.size(); ++i) {
            file << parameters[i];
            if (i < parameters.size() - 1)
                file << "\t";  // Tab-separated columns
        }
        file << "\n";

        // Write values to the file
        for (size_t i = 0; i < values.size(); ++i) {
            for (size_t j = 0; j < values[i].size(); ++j) {
                file << values[i][j];
                if (j < values[i].size() - 1)
                    file << "\t";
            }
            file << "\n";
        }

        file.close(); // Close file
    }
//generates combination excluding permutations of single particle states (quantum number triplets)
std::vector<std::array<int, 6>>* generateCombinations(int nmin, int umin, int mmin, int nmax, int umax, int mmax) { 
    auto combinations = new std::vector<std::array<int, 6>>;

    for (int n1 = 0; n1 <= nmax; ++n1) {
        for (int u1 = 0; u1 <= umax; ++u1) {
            for (int m1 = 0; m1 <= mmax; ++m1) {
                for (int n2 = 0; n2 <= nmax; ++n2) {
                    for (int u2 = 0; u2 <= umax; ++u2) {
                        for (int m2 = 0; m2 <= mmax; ++m2) {
                            std::array<int,6> target_array = {n2,u2,m2,n1,u1,m1};
                            if (std::find(combinations->begin(), combinations->end(), target_array) != combinations->end()) {
                                    continue;
                                }//reduces output array by half (doesn't allow repeating permutations of multiindeces)

                            std::array<int,6> new_element = {n1, u1, m1, n2, u2, m2};
                            
                            combinations->push_back(new_element);
                        }
                    }
                }
            }
        }
    }

    return combinations;
}
void print_matrix(double** matrix, unsigned int size){
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << matrix[i][j] << ", ";
        }
        std::cout << std::endl;
    }
}

void write_to_csv(const std::string &filename, 
                  const std::vector<double> &x, 
                  const std::vector<double> &y,
                  int nmax, int umax, int mmax) {
    // Open the file
    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::invalid_argument("Unable to open the file: " + filename);
        return;
    }
    file << "nmax: " << nmax << ", umax: " << umax << ", mmax: " << mmax << "\n";
    file << "X,Y\n";
    if (x.size() != y.size()) {
        throw std::invalid_argument("The two arrays have different lengths.");
        file.close();
        return;
    }

    // Write data
    for (size_t i = 0; i < x.size(); ++i) {
        file << x[i] << "," << y[i] << "\n";
    }

    // Close the file
    file.close();
}