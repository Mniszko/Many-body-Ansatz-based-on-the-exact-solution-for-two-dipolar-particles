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
#include <gsl/gsl_integration.h>
#include <fstream>

double h=1.;

double fast_power(double a, int b){
    double result = 1;
    for(int i=0; i<b ; ++i){
        result *= a;
    }
    return result;
}
// Simple RAII wrapper
class IntegrationWorkspace {
  gsl_integration_workspace *wsp;

public:
  IntegrationWorkspace(const size_t n = 1000)
      : wsp(gsl_integration_workspace_alloc(n)) {}
  ~IntegrationWorkspace() { gsl_integration_workspace_free(wsp); }

  operator gsl_integration_workspace *() { return wsp; }
};

// Build gsl_function from lambda
template <typename F> class gsl_function_pp : public gsl_function {
  const F func;
  static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp *>(params)->func(x);
  }

public:
  gsl_function_pp(const F &f) : func(f) {
    function = &gsl_function_pp::invoke; // inherited from gsl_function
    params = this;                       // inherited from gsl_function
  }
  operator gsl_function *() { return this; }
};

// Helper function for template construction
template <typename F> gsl_function_pp<F> make_gsl_function(const F &func) {
  return gsl_function_pp<F>(func);
}


/*dispersion test*/
/*
int main(){
    int nmin=0;
    int umin=0;
    int mmin=0;
    int nmax=2;
    int umax=1;
    int mmax=1;

    std::vector<std::array<double,2>> norms;

    double epsabs = 1e-8;
    double epsrel = 1e-6;
    size_t limit = 100;

    char* argv[2];
    double length = 0.01;

    CompleteIntegral integral(limit,epsabs,epsrel);
    std::vector<std::array<int,6>> combination = *generateCombinations(nmin,umin,mmin,nmax,umax,mmax);
    unsigned int number_of_functions = combination.size();
    main_loop(1, argv, 1, 1, length, nmax, umax, mmax, true);
    double* array = readLastColumn("./vector_outputs/lowest_wavefunction"+std::to_string(length)+".csv", number_of_functions);
    Eigen::MatrixXd changerr(number_of_functions,number_of_functions);
    Eigen::MatrixXd changerr2(number_of_functions,number_of_functions);
    integral.change_length(1);
    for (int i=0 ; i < number_of_functions ; ++i){
        for (int j=0 ; j <= i ; ++j){
            std::array<int,6> l1l2 = combination.at(i);
            std::array<int,6> l3l4 = combination.at(j);

            double c1=array[i];
            double c2=array[j];

            //double c1=array[i];
            //double c2=array[j];

            integral.changeMultiIndex(l1l2[0], l1l2[1], l1l2[2], 1);
            integral.changeMultiIndex(l1l2[3], l1l2[4], l1l2[5], 2);
            integral.changeMultiIndex(l3l4[0], l3l4[1], l3l4[2], 3);
            integral.changeMultiIndex(l3l4[3], l3l4[4], l3l4[5], 4);

            integral.change_length(length);


            if (std::abs(c1*c2)<1e-8){ // there are nearly nonexisting coefficients that just bloat the calculations
                changerr(i,j) = 0;
                changerr2(j,i) = 0;
                changerr(j,i) = 0;
                changerr2(i,j) = 0;
            }


            double normi = 1 / std::sqrt(integral.integrate_norm_external(l1l2[0], l1l2[1], l1l2[2]));
            double normj = 1 / std::sqrt(integral.integrate_norm_external(l1l2[3], l1l2[4], l1l2[5]));
            double normk = 1 / std::sqrt(integral.integrate_norm_external(l3l4[0], l3l4[1], l3l4[2]));
            double norml = 1 / std::sqrt(integral.integrate_norm_external(l3l4[3], l3l4[4], l3l4[5]));

            double r1_val = integral.integrate_r1(normi,normj,normk,norml);
            double r1_squared_val = integral.integrate_r1_squared(normi,normj,normk,norml);
            double r1r2_val = integral.integrate_r1r2(normi,normj,normk,norml);

            // reason for using this code would be additionally getting only absolute values of c values.
            //r1_val = std::abs(r1_val);
            //r1_squared_val = std::abs(r1_squared_val);
            //r1r2_val = std::abs(r1r2_val);

            //explicit to be understandable
            // symnorm equals to 1/2 or 1/sqrt(2) for each wavefunction
            double r_elem = c1*c2*(r1_val+r1_val)*0.5;
            double r_squared_elem = c1*c2*(2*r1_squared_val+2*r1r2_val)*0.25;
            r_squared_elem=c1*c2*r1_squared_val;
            
            if (i==0 && j==0){
                std::cout << c1 << ' ' << c2 << std::endl;
                std::cout << "squared: " << r1_squared_val << std::endl;
                std::cout << "single: " << r1_val << std::endl;
                std::cout << l1l2[0] << l1l2[1] << l1l2[2] << std::endl;
                std::cout << l1l2[3] << l1l2[4] << l1l2[5] << std::endl;
                std::cout << l3l4[0] << l3l4[1] << l3l4[2] << std::endl;
                std::cout << l3l4[3] << l3l4[4] << l3l4[5] << std::endl;
            }
            changerr(i,j) = r_elem;
            changerr2(j,i) = r_squared_elem;
            changerr(j,i) = r_elem;
            changerr2(i,j) = r_squared_elem;
        }
    }
    Eigen::MatrixXd changer =changerr;

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

    file.close();

    changer =changerr2;

    file.open("./testy/sample_matrix2.csv");
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

    double rsum = changerr.sum();
    double r2sum = changerr2.sum();

    //rsum=r_matrix(0,0);
    //r2sum = r2_matrix(0,0);
    std::cout << "r2sum = " <<r2sum << std::endl;
    std::cout << "rsum = " << rsum << std::endl;
    double radial_dispersion = sqrt(r2sum-rsum*rsum);
    std::cout << "dyspersja = " << radial_dispersion << std::endl;
    std::cout << "dyspersja w potędze = " << r2sum-rsum*rsum << std::endl;
    return 0;
}
*/

int main(){
    double epsabs = 1e-8;
    double epsrel = 1e-6;
    size_t limit = 100;

    int key = 6;

    //CompleteIntegral integral(limit,epsabs,epsrel);

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

    //double norm1 = integral.integrate_norm_external(n1, u1, m1);
    //double norm2 = integral.integrate_norm_external(n2, u2, m2);
    //double norm3 = integral.integrate_norm_external(n3, u3, m3);
    //double norm4 = integral.integrate_norm_external(n4, u4, m4);

    //integral.changeMultiIndex(n1, u1, m1, 1);
    //integral.changeMultiIndex(n2, u2, m2, 2);
    //integral.changeMultiIndex(n3, u3, m3, 3);
    //integral.changeMultiIndex(n4, u4, m4, 4);

    double length = 1;
    double rev_length = 1/length;
    double result, abserr, outer_result, outer_abserr, middle_result,
      middle_abserr, inner_result, inner_abserr;

    IntegrationWorkspace wsp1(limit);
    IntegrationWorkspace wsp11(limit);
    IntegrationWorkspace wsp2(limit);
    IntegrationWorkspace wsp22(limit);
    IntegrationWorkspace wsp3(limit);
    IntegrationWorkspace wsp33(limit);
    IntegrationWorkspace wsp4(limit);
    IntegrationWorkspace wsp44(limit);

    double epsilon = 2e-3;

    int smthng = 0;
    int highsmthng = 0;
    double templimit = 15;
    double pts[3];
    pts[0] = 0, pts[2]=20;
    auto full = make_gsl_function([&](double z2) {
        if (highsmthng==limit){
            std::cout << "made progress: " << highsmthng << "\taround " << smthng << " to go" << std::endl;
            limit += limit;
        }
        smthng =0;
        ++highsmthng;
        auto outer = make_gsl_function([&](double z1) {
            ++smthng;
        auto middle = make_gsl_function([&](double rho2) {
            auto inner = make_gsl_function([&](double rho1) {
            // for now integration only for two particles in state |j,i>=|l,k>=|0(3),0(3)>
                        
            double sillyAndGoofy = (rho2-rho1)*(rho2-rho1)+(z2-z1)*(z2-z1) + epsilon;
            
            //this may produce unpredictable behavior
            //to somehow justify this - we will also calculate delta potential ensuring hard-core behavior. I don't know if this will be sufficient, but I am not familiar with methods of integrating such singularities in so many dimensions.
            
            
            double exponential = std::exp(-rho1*rho1-rho2*rho2);
            double cosSqTheta = (z2-z1)*(z2-z1)/(sillyAndGoofy);
            double evil = fast_power(std::sqrt(sillyAndGoofy), 5);
            pts[1] = rho2; 
            return exponential*(1-3*cosSqTheta)/evil;

            });
            // qags because we have singularity on 0
            double addon_result, addon_abserr;
            gsl_integration_qagp(inner, pts, 3, epsabs, epsrel, limit, wsp1,
                                &inner_result, &inner_abserr);
            //gsl_integration_qagiu(inner, epsilon, epsabs, epsrel, limit, wsp11, &addon_result, &addon_abserr);
            return inner_result + addon_result;
        });
        double addon_result, addon_abserr;
        gsl_integration_qag(middle, 0, 20, epsabs, epsrel, limit, key, wsp2,
                            &middle_result, &middle_abserr);
        //gsl_integration_qagiu(middle, epsilon, epsabs, epsrel, limit, wsp11,&addon_result, &addon_abserr);
        return middle_result + addon_result;
        });
        double addon_result, addon_abserr;
        gsl_integration_qag(outer, 0, length, epsabs, epsrel, limit, key,
                            wsp3, &outer_result, &outer_abserr);
        //gsl_integration_qagiu(outer, epsilon, epsabs, epsrel, limit, wsp11,&addon_result, &addon_abserr);
        return outer_result + addon_abserr;
    });
    double addon_result, addon_abserr;
    gsl_integration_qag(full, 0, length, epsabs, epsrel, limit, key, wsp3,
                        &result, &abserr);
    //gsl_integration_qagiu(full, epsilon, epsabs, epsrel, limit, wsp11,&addon_result, &addon_abserr);
        
    std::cout << "result: " << result << std::endl;
     
    return 0;
}
/*
int main(){
    double result, abserr;
    double epsabs = 1e-8;
    double epsrel = 1e-6;
    size_t limit = 100;

    int key = 6;

    IntegrationWorkspace wsp1(limit);

    double Z2MinusZ1 = 0;
    double epsilon = 1e-8;
    //test sampling integration to test out ways of getting away with singularities
    //and using different integration methods
    auto integrand = make_gsl_function([&](double rho) {

        double exponential = std::exp(-rho*rho);
        double evil = sqrt(rho*rho+Z2MinusZ1)+epsilon;
        double potentialTerm = (1-3*Z2MinusZ1/evil)/fast_power(evil, 5); 
        double finval = exponential*potentialTerm;



        return finval;
    });

    std::cout << "qag method" << std::endl;
    gsl_integration_qag(integrand, -1, 20, epsabs, epsrel, limit, key, wsp1,&result, &abserr);
    std::cout << "result: "<< result << std::endl;
    std::cout << "error: " << abserr << std::endl;

    std::cout << "qags method" << std::endl;
    gsl_integration_qags(integrand, -1, 20, epsabs, epsrel, limit, wsp1,&result, &abserr);
    std::cout << "result: "<< result << std::endl;
    std::cout << "error: " << abserr << std::endl;

    std::cout << "qagiu method" << std::endl;
    gsl_integration_qagiu(integrand, 0, epsabs, epsrel, limit, wsp1,&result, &abserr);
    std::cout << "result: "<< result << std::endl;
    std::cout << "error: " << abserr << std::endl;
    return 0;
}
*/
/*

we're using limited range for rho due to disperssion << 10
gaussian interpolation is set to most robust (key=6)

without singularities qags and qag methods are the most fit, but all of the above achieve errors of simmilar order.

integrating only gaussian works best with qags surprisingly

qags and qagiu can't deal with very small values in devisor

THEOREM: singularities that appear in this equation can not be integrated numerically with nested integrals.
REASONING: integrating separately through two symmetrically defined variables will result in one specific integral giving infinite result while the other one should restrict the value going to infinity (see plots of specific). With this reason one can not simply integrate through both variables one after another, instead all singularities that cancel one another should be integrated at once (using qags method). Either pairs of integrals should have volume forms given by drdz or dzdz and drdr.

other approaches I will try include:
-integrating with positive constant epsilon added to the devisor
[this method fails as results differ strongly due to epsilon values or takes too much time to be calculated for lower values. qags, qagp and qag all the same can't manage with integrating that abomination]
-integrating with "hard core" approximation achieved by equating integrand to 0 for low values of devisor
[this bares no useful result. Both qag and qags algorithms fail to integrate that function giving off errors]

[thus any hard defined integration (e.i. trapezoidal, rectangular) and even recursive interpolation will probably fail the task]

plan for now:
- develop better reasoning behind mentioned theorem
- try to  actually do something with it

if I won't succeed I will just replicate Michał's results for bigger matrices

*/
