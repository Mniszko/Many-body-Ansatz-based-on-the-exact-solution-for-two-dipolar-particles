#include "variables.h"

#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <array>
#include <iostream>
#include <boost/program_options.hpp>
#include <stdexcept>
#include <dirent.h>
#include <unistd.h>

const double pi = 3.141592653589793;

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
std::vector<std::array<int, 6>>* generateCombinations(int nmin, int umin, int mmin, int nmax, int umax, int mmax, bool positive, double cutoff) { 
    auto combinations = new std::vector<std::array<int, 6>>;
    if (nmax!=-1){
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
                                if (positive){
                                    std::array<int,6> new_element = {n1, u1, m1, n2, u2, m2};
                                    combinations->push_back(new_element);
                                } else if (n1-n2==0 && m1-m2 == 0) {
                                    
                                    std::array<int,6> new_element = {n1, u1, m1, -n2, u2, -m2};
                                    combinations->push_back(new_element);
                                } else {
                                    continue;
                                }

                            }
                        }
                    }
                }
            }
        }
    } else if (nmax==-1 && umax==-1 && mmax==-1){//activates energy cutoff is parsed as argument of main (cutoff is stored as mmax)
        double cutoff_energy = cutoff;
        std::vector<std::array<int, 3>> triplets;
        double length = 7;
        for (int n = 15; n >= 0; n--) {
            for (int u = 15; u >= 0; u--) {
                for (int m = u; m >= 0; m--) {
                    double E = 2*(pi * pi) / length/length * n * n + 2 * u + m + 1;

                    if (E < cutoff_energy) {
                        triplets.push_back({n, u, m});
                    }
                }
            }
        }
        for (size_t i = 0; i < triplets.size(); ++i) {
            for (size_t j = 0; j < triplets.size(); ++j) {
                std::array<int,6> new_element = {triplets[i][0],triplets[i][1],triplets[i][2],-triplets[j][0],triplets[j][1],-triplets[j][2]};
                
                //if (std::find(combinations->begin(), combinations->end(), new_element) != combinations->end()) {
                //        continue;
                //    }
                if (triplets[i][0] == triplets[j][0]){
                    combinations->push_back(new_element);
                }
            }
        }
    } else {
        std::cerr << "error occured while generating vectors" <<std::endl;
    }
    if (combinations->size() == 0){
        std::cerr << "matrix has no elements" << std::endl;
        exit(1);
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


void write_vectors_to_csv(const std::vector<std::array<int, 6>>& vector_of_vectors, 
                double* stdArray, 
                size_t arraySize, 
                const std::string& filename) {

    // Check if sizes match; if not, something is wrong
    if(vector_of_vectors.size() != arraySize) {
        std::cerr << "Mismatch in sizes of vector_of_vectors and stdArray!" << std::endl;
        return;
    }

    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    file << "n1,u1,m1,n2,u2,m2,constant\n";
    for (size_t i = 0; i < vector_of_vectors.size(); ++i) {
        for (int j = 0; j < 6; ++j) {
            file << vector_of_vectors[i][j];
            if (j != 5) { // Don't put a comma after the last value in the array
                file << ",";
            }
        }
        // After writing the contents of vector_of_vectors, write the corresponding value from stdArray
        file << "," << stdArray[i];
        file << "\n"; // New line after each row
    }

    file.close();
}

double* readLastColumn(const std::string& filename, unsigned int& size) {
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return nullptr;
    }

    std::string line;
    
    // Skip header
    std::getline(file, line);

    std::vector<double> values;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        double lastValue = 0.0;

        while (std::getline(iss, token, ',')) {
            try {
                lastValue = std::stod(token);
            } catch (const std::invalid_argument&) {
                continue; // handle invalid conversion, move to the next token
            }
        }
        values.push_back(lastValue);
    }

    size = values.size();
    double* result = new double[size];
    
    for (size_t i = 0; i < size; ++i) {
        result[i] = values[i];
    }

    return result;
}

void deleteFilesWithPrefix() {
    const std::string directory = "./vector_outputs/";
    const std::string prefix = "lowest_wavefunction";
    DIR* dir;
    struct dirent* ent;
    int numbOfFiles = 0;
    if ((dir = opendir(directory.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (strncmp(ent->d_name, prefix.c_str(), prefix.size()) == 0) {
                std::string filePath = directory + "/" + ent->d_name;
                if (unlink(filePath.c_str()) == 0) {
                    ++numbOfFiles;
                } else {
                    std::cerr << "Error deleting: " << filePath << std::endl;
                }
            }
        }
        closedir(dir);
    } else {
        std::cerr << "Could not open directory: " << directory << std::endl;
    }
    std::cout << "Deleted " << std::to_string(numbOfFiles) << " files" << std::endl;
}

