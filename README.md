# >>software nickname<<

## Introduction

This software is intended to calculate Many-body Ansatz based on the exact solution for two dipolar particles in 3-dimensional trap (cylindrically symmetric harmonic oscillator with longitudinal dimension being a box with periodic boundary conditions) using numerical methods on multiple quantum states.

## Structure

- Repository directory contains main file for running software (g-plot.cpp), main file for testing software (test.cpp), makefile building executables *g-plot* or *test* and (after successfully calculating outputs) tables with extension *.csv* that store values that can be plotted with *plot.py*
- vector_outputs stores vectors corresponding with lowest energy state for 2 particle system encoded in a base of two single particle eigenfunctions $\lvert (n,u,m)_1,(n,u,m)_2\rangle$ where $n$, $u$ and $m$ are longitudinal quantum number, radial quantum number and angular quantum number respectively. Last row contains normalized coefficients of lowest energy vector in noninteracting system basis.
- src directory contains Eigen directory (explained in next section), headers and definitions of classes and functions.
    - *integral.h* contains declaration of CompleteIntegral class storing quantum numbers to make calculating integrals for multiple matrix elements intuitive
    - *flag_parser.h* contains declaration of class taking care of flaggs being passed to executable (currently its purpose is to help with debugging, it will change in future)
    - *main_functions.h* contains declaration of functions that are directly called by *main* function - so all necesary loops over multiple matrix elements and loops calling previous ones with different parameters to generate output
    - *variables.h* contains helper functions for writing data into files, reading data, generating encoding for eigenvectors, etc.

## Dependencies (and installation)
For matrix manipulation [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library is used. To properly install Eigen on linux, one can use following snippet from :
```bash
wget -O Eigen.zip https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
unzip Eigen.zip
mv ./eigen-3.4.0/Eigen ./src/
rm Eigen.zip
rm -r ./eigen-3.4.0/
```
or download stable release (for example 3.4.0) from [Eigen site](https://eigen.tuxfamily.org/index.php?title=Main_Page) and extract Eigen directory to ./src/ directory in repository directory.

Gnu Scientific Library ([GSL](https://www.gnu.org/software/gsl/doc/html/integration.html)) for C++ takes care of integration. It can be installed via standard package managers on linux or using [this](https://www.gnu.org/software/gsl/extras/native_win_builds.html) or [this](https://solarianprogrammer.com/2020/01/26/getting-started-gsl-gnu-scientific-library-windows-macos-linux/) tutorial for Windows

## Integration techniques
for integration without singularities following functions are used:
- over infinite intervals from 0 to $\infty$: **gsl_integration_qagiu** which mapps function to finite interval and applies Gauss-Kronrod 21-point integration rule adaptively until errors for integral are below set threshold (defined with integral class defined in the file *integral.h*)
- over finite integrals: **gsl_integration_qag** with predefined key equal to 6 (so 61 point Gauss-Kronrod rules), which can be changed with change_key function in class CompleteIntegral inside file *integral.h*. Key is set to maximum due to very low accuracy for lower keys (especially key equal to 1), but lowering it may speed up the code with maintained high accuracy.

### Implementation of integration techniques:
Integration constants are defined in CompleteIntegral class constructor.

Integration workspaces storing necessary data for integration process are constructed only after calling specific functions inside CompleteIntegral class.

Higher order integration is achieved by nesting integrals inside one another - for this purpose helper functions and classes are defined before defining specific integration functions. This implementation is based on [this](https://stackoverflow.com/a/43636411/12955940) answer to forum post.

## Usage

To run specific loop one has to uncomment it (and prefferably comment out the rest) in *g-plot.cpp*, change maximum quantum numbers in the same file (for example nmax=2 means that generated states include $n\in\{0,1,2\}$), build the executable *g-plot* with command
```bash
make
```
and run it.

Plotting can be done with preferred software by using *.csv* outputs or with python using *plot.py* and commenting out proper function at the end of the file.

Currently to change initial and maximum values for length $L$ and variable $g_\delta$ it is required to do so inside specific function in *main_funcftions.cpp*. For example to calculate 100 lowest energy points between $L=1$ and $L=10$ it is necessary to set:
```python
(...)
    double lengthmax = 10;
    const unsigned int numofsteps = 100;
    (...)
    for (double length=1; length < lengthmax ; length+=lengthmax/numofsteps){
(...)
```
inside *length_loop* function.

## Limitations

For very low and very high longitudinal lengths $L$ predictions are somewhat inaccurate - in both of these cases orthonormal interacting 2-boson wavefunctions are made of multiple orthonormal noninteracting 2-boson wavefunctions with high coefficients. This software is capable of calculating simultaneously only limited (up to around 8 for each particle with reasonable time of computation) number of states, thus its predictions should be useful for medium lengths $L$ - (probably) around 1-10 units.

As for time - every iteration of disperion loop for 8 one particle states takes around 10 seconds - it can be made significantly faster with additional rules for integration (most integrals are being repeated).  Calculation time for energy matrix scales up very fast, so one has to be very careful with how many states are taken into account.

Currently integrating only interactions in cylindrical system is possible due to multiple shortcuts being utilized to make the code faster. CompleteIntegral class can be easily expanded, but currently doesn't support general way of integrating interactions or dispersion between particles in custom traps.

## To do

- [x] developing universal code structure, dividing program into specialized files
- [x] energy matrix for hard-core delta interactions and its diagonalization
- [x] exporting results in .csv format and plotting them
- [x] saving full vectors and matrices, making testfile for troubleshooting purposes
- [x] extracting data for lowest energy state in function of $g_\delta$ and $L$
- [x] extracting data for highest coefficient in function of $L$
- [x] extracting data for radial dispersion in terms of geometric center (observable $\hat{\rho}_1$) and center of mass (observable $(\hat{\rho}_1+\hat{\rho}_2)/2$)
- [ ] energy matrix for dipole interactions and its diagonalization with 4'th rank integrals
- [ ] energy matrix for dipole interactions and its diagonalization with 4'th rank integrals with shorthands for specific states
- [ ] energy matrix for dipole interactions and its diagonalization with 3'th rank integrals (maybe with shorthands for specific states)
- [ ] extracting lowest energy and radial dispersion for combined delta and dipole interactions, in their respectable functions
- [ ] writing combined loop for running all calculations at the same time (it will make the final code two- to three-fold faster)
- [ ] use flag parser to update maximum quantum numbers and select loops by parsed flags

## References
- efficient integration and code simplification using GSL (helper functions and classes for writing integrals) - [here](https://stackoverflow.com/a/43636411/12955940)