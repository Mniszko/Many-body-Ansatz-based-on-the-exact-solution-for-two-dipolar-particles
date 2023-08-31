#include "integral.h"

#include "../src/Eigen/Dense"
#include "../src/Eigen/Eigenvalues"
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>


const long double hbar = 1; /*1.0545718e-34;*/
const double pi = 3.14159265358979323846;
const double sqrtPi = 1.7724538509055159;

//const double length = 1; /*długość małego walca == L*/
//const double rev_length = 1/length;
const double radius = 1; /*promień małego walca == sigma*/
const double mass = 1;

const int NumberOfSteps = 8e2;
const double dstep=1./double(NumberOfSteps);
const double h = 1;

// Simple RAII wrapper 
class IntegrationWorkspace {
    gsl_integration_workspace * wsp;

    public:
    IntegrationWorkspace(const size_t n=1000):
        wsp(gsl_integration_workspace_alloc(n)) {}
    ~IntegrationWorkspace() { gsl_integration_workspace_free(wsp); }

    operator gsl_integration_workspace*() { return wsp; }
};

// Build gsl_function from lambda
template <typename F>
class gsl_function_pp: public gsl_function {
    const F func;
    static double invoke(double x, void *params) {
        return static_cast<gsl_function_pp*>(params)->func(x);
    }
    public:
        gsl_function_pp(const F& f) : func(f) {
        function = &gsl_function_pp::invoke; //inherited from gsl_function
        params   = this;                     //inherited from gsl_function
        }
        operator gsl_function*(){return this;}
};

// Helper function for template construction
template <typename F>
gsl_function_pp<F> make_gsl_function(const F& func) {
    return gsl_function_pp<F>(func);
}

void CompleteIntegral::check_if_valid_key(int k){
    if(k != 1 && k != 2 && k != 3 && k != 4 && k != 5 && k != 6){
                throw std::runtime_error("Wrong value for gaussian interpolation key");
            }
}
double CompleteIntegral::fast_power(double a, int b){
    double result = 1;
    for(int i=0; i<b ; ++i){
        result *= a;
    }
    return result;
}

void CompleteIntegral::print_result(){
            std::cout << "result:\t" << this->result << "\nerror:\t" << this->error << std::endl;
        }

void CompleteIntegral::changeMultiIndex(int n, int u, int m, int changer){
    if (changer==1){
        this->l1 = {n, u, m};
    } else if (changer==2){
        this->l2 = {n, u, m};
    } else if (changer==3){
        this->l3 = {n, u, m};
    } else if (changer==4){
        this->l4 = {n, u, m};
    }
}
void CompleteIntegral::changeLimit(int limit){
    this->limit = limit;
}
double CompleteIntegral::radial_function(double r, int n, int u, int m){
    double psi = std::exp(-h*r*r*0.5)*this->fast_power(r,m)*boost::math::hypergeometric_1F1(-u, m + 1, r*r * h);
    return psi;
}

//czy to jest w ogóle potrzebne?
/*
void CompleteIntegral::integrate_norm() {
    double result1, error1, result2, error2;
    IntegrationWorkspace wsp1(this->limit);

    gsl_function F = make_gsl_function( [&](double rho) {
            // inside of a function
            double exponent = -h * rho * rho;
            double r_power = std::pow(rho, (double)2 * this->m1);
            double hypergeometric1F1 = boost::math::hypergeometric_1F1(-this->u1, 1 + this->m1, h * rho * rho);

            return exp(exponent) * r_power * hypergeometric1F1 * hypergeometric1F1; 
        } );
    gsl_integration_qagiu(&F, 0, this->epsabs, this->epsrel, this->limit, wsp1, &result1, &error1);

    if (this->n1!=this->n2 || this->m1!=this->m2 || this->u1!=this->u2){
        IntegrationWorkspace wsp2(this->limit);
        gsl_function F = make_gsl_function( [&](double rho) {
            // inside of a function
            double exponent = -h * rho * rho;
            double r_power = pow(rho, 2 * this->m2);
            double hypergeometric1F1 = boost::math::hypergeometric_1F1(-this->u2, 1 + this->m2, h * rho * rho);

            return exp(exponent) * r_power * hypergeometric1F1 * hypergeometric1F1; 
        } );
        gsl_integration_qagiu(&F, 0, this->epsabs, this->epsrel, this->limit, wsp2, &result2, &error2);
    } else {
        result2 = result1;
        error2 = error1;
    }
    this->norm1 = result1*2*pi;
    this->norm2 = result2*2*pi;
}
*/
double CompleteIntegral::integrate_norm_external(int n,int u,int m) {
    double result, error;
    IntegrationWorkspace wsp(this->limit);

    gsl_function F = make_gsl_function( [&](double rho) {
            // inside of a function
            double exponent = -h * rho * rho;
            double r_power = std::pow(rho, (double)2 * m);
            double hypergeometric1F1 = boost::math::hypergeometric_1F1(-u, 1 + m, h * rho * rho);

            return exp(exponent) * r_power * hypergeometric1F1 * hypergeometric1F1; 
        } );
    gsl_integration_qagiu(&F, 0, this->epsabs, this->epsrel, this->limit, wsp, &result, &error);
    return result;
}


double CompleteIntegral::fast_add_over_harmonic(unsigned int number_of_first_function, unsigned int number_of_second_function){
    if (number_of_first_function!=number_of_second_function){
        return 0;
    }
    int n1 = this->l1[0];
    int u1 = this->l1[1];
    int m1 = this->l1[2];
    int n2 = this->l2[0];
    int u2 = this->l2[1];
    int m2 = this->l2[2];
    double result = (n1*n1+n2*n2)*0.5*pi*pi*this->rev_length*this->rev_length + 2*u1+2*u2+m1+m2+2;
    return result;
}

double CompleteIntegral::integrate_over_delta(double g_del){
    double result, abserr;
    auto l = this->l1;
    int n1 = l[0];
    int u1 = l[1];
    int m1 = l[2];
    l = this->l2;
    int n2 = l[0];
    int u2 = l[1];
    int m2 = l[2];
    l = this->l3;
    int n3 = l[0];
    int u3 = l[1];
    int m3 = l[2];
    l = this->l4;
    int n4 = l[0];
    int u4 = l[1];
    int m4 = l[2];
    if (m3*m4!=m1*m2 || n3*n4!=n1*n2){
        return 0;
    }
    IntegrationWorkspace wsp(this->limit);
    auto inner = make_gsl_function( [&](double rho) {
                    // inside of a function
                    /*
                    double psi1 = std::exp(-h*rho*rho*0.5)*this->fast_power(rho,m1)*boost::math::hypergeometric_1F1(-u1, m1 + 1, rho*rho * h);
                    double psi2 = std::exp(-h*rho*rho*0.5)*this->fast_power(rho,m2)*boost::math::hypergeometric_1F1(-u2, m2 + 1, rho*rho * h);
                    double psi3 = std::exp(-h*rho*rho*0.5)*this->fast_power(rho,m3)*boost::math::hypergeometric_1F1(-u3, m3 + 1, rho*rho * h);
                    double psi4 = std::exp(-h*rho*rho*0.5)*this->fast_power(rho,m4)*boost::math::hypergeometric_1F1(-u4, m4 + 1, rho*rho * h);
                    */
                    double psi1 = radial_function(rho, n1, u1, m1);
                    double psi2 = radial_function(rho, n2, u2, m2);
                    double psi3 = radial_function(rho, n3, u3, m3);
                    double psi4 = radial_function(rho, n4, u4, m4);
                    return psi1*psi2*psi3*psi4;
                } );
    gsl_integration_qagiu(inner, 0, epsabs, epsrel, limit, wsp,
                                    &result, &abserr);
    return result*g_del;
}

double CompleteIntegral::dipole_integral_function(std::array<int,3> l1, std::array<int,3> l2, std::array<int,3> l3, std::array<int,3> l4){
    int n1 = l1[0];
    int u1 = l1[1];
    int m1 = l1[2];

    int n2 = l2[0];
    int u2 = l2[1];
    int m2 = l2[2];

    int n3 = l3[0];
    int u3 = l3[1];
    int m3 = l3[2];

    int n4 = l4[0];
    int u4 = l4[1];
    int m4 = l4[2];
    double result, abserr, outer_result, outer_abserr, middle_result, middle_abserr, inner_result, inner_abserr;

    IntegrationWorkspace wsp1(this->limit);
    IntegrationWorkspace wsp11(this->limit);
    IntegrationWorkspace wsp2(this->limit);
    IntegrationWorkspace wsp22(this->limit);
    IntegrationWorkspace wsp3(this->limit);
    IntegrationWorkspace wsp33(this->limit);
    IntegrationWorkspace wsp4(this->limit);
    IntegrationWorkspace wsp44(this->limit);

    auto full= make_gsl_function([&](double z2){
        auto outer = make_gsl_function([&](double z1){
            auto middle = make_gsl_function( [&](double rho2) {
                auto inner = make_gsl_function( [&](double rho1) {
                    //I'm not sure if only cosines will do, although I have gut feeling they will
                    double psi1 = std::exp(-h*rho1*rho1*0.5)*this->fast_power(rho1,m1)*boost::math::hypergeometric_1F1(-u1, m1 + 1, rho1*rho1 * h);
                    double Z1 = std::cos(2*pi*this->rev_length*n1*z1);
                    double psi2 = std::exp(-h*rho2*rho2*0.5)*this->fast_power(rho2,m2)*boost::math::hypergeometric_1F1(-u2, m2 + 1, rho2*rho2 * h);
                    double Z2 = std::cos(2*pi*this->rev_length*n2*z2);
                    double psi3 = std::exp(-h*rho1*rho1*0.5)*this->fast_power(rho1,m3)*boost::math::hypergeometric_1F1(-u3, m3 + 1, rho1*rho1 * h);
                    double Z3 = std::cos(2*pi*this->rev_length*n3*z1);
                    double psi4 = std::exp(-h*rho2*rho2*0.5)*this->fast_power(rho2,m4)*boost::math::hypergeometric_1F1(-u4, m4 + 1, rho2*rho2 * h);
                    double Z4 = std::cos(2*pi*this->rev_length*n4*z2);
                    
                    double distance = sqrt((rho2-rho1)*(rho2-rho1) + (z2-z1)*(z2-z1));
                    double dipole = ((rho2-rho1)*(rho2-rho1) - (z2-z1)*(z2-z1))/(distance*distance*distance*distance*distance);

                    return psi1*psi2*psi3*psi4*dipole*Z1*Z2*Z3*Z4;
                });
                //qags because we have singularity on 0
                double epsilon = 1e-7;
                double addon_result, addon_abserr;
                gsl_integration_qags(inner, 0, epsilon, epsabs, epsrel, limit, wsp1,
                                        &inner_result, &inner_abserr);
                gsl_integration_qagiu(inner, epsilon, epsabs, epsrel, limit, wsp11, &addon_result, &addon_abserr);
                return inner_result+addon_result;
            });
            double epsilon = 1e-7;
            double addon_result, addon_abserr;
            gsl_integration_qags(middle, 0, epsilon, epsabs, epsrel, limit, wsp2,&middle_result, &middle_abserr);
            gsl_integration_qagiu(middle, epsilon, epsabs, epsrel, limit, wsp11, &addon_result, &addon_abserr);
            return middle_result+addon_result;
        });
        double epsilon = 1e-7;
        double addon_result, addon_abserr;
        int key=this->key;
        gsl_integration_qag(outer, 0, this->length, epsabs, epsrel, limit, key, wsp3, &outer_result, &outer_abserr);
            gsl_integration_qagiu(outer, epsilon, epsabs, epsrel, limit, wsp11, &addon_result, &addon_abserr);
        return outer_result+addon_abserr;
    });
    double epsilon = 1e-7;
    double addon_result, addon_abserr;
    int key=this->key;
    gsl_integration_qag(full, 0, this->length, epsabs, epsrel, limit, key, wsp3, &result, &abserr);
            gsl_integration_qagiu(full, epsilon, epsabs, epsrel, limit, wsp11, &addon_result, &addon_abserr);

    return result+addon_abserr;
}

// 4-fold integral, the slowest one
double CompleteIntegral::integrate_over_dipole(double g_dip){
    std::array<int,3> l1 = this->l1, l2 = this->l2, l3 = this->l3, l4 = this->l4;
    int u1 = l1[1];
    int m1 = l1[2];

    int u2 = l2[1];
    int m2 = l2[2];

    int u3 = l3[1];
    int m3 = l3[2];

    int u4 = l4[1];
    int m4 = l4[2];

    if ((m1!=m3 && m2!=m4)||((m2!=m3 && m1!=m4))){
        return 0;
    } // later I'll add extra if statement for second order integrals - if I find proper one
    // besides there may be a quicker solution for specific quantum numbers, it would be great to find one
    else if (m3==m1 && m4==m2 && m3!=m4){
        double result = 0;
        result += this->dipole_integral_function(l1, l2, l4, l3);
        result += this->dipole_integral_function(l2, l1, l3, l4);
        return result;
    }
    else if (m4==m1 && m3==m2 && m3!=m4){
        double result = 0;
        result += this->dipole_integral_function(l2, l1, l4, l3);
        result += this->dipole_integral_function(l1, l2, l3, l4);
        return result;
        }
    else { //only m1=m2=m3=m4
        double result = 0;
        result += this->dipole_integral_function(l1, l2, l4, l3);
        result += this->dipole_integral_function(l2, l1, l3, l4);
        result += this->dipole_integral_function(l2, l1, l4, l3);
        result += this->dipole_integral_function(l1, l2, l3, l4);
        return result;
    }


    return result*g_dip;
}

double CompleteIntegral::integrate_delta_dispersion_helper_r(double c1, double c2, bool centerOfMass){
    double result_inner, abserr_inner, result, abserr;
    double epsabs=this->epsabs, epsrel=this->epsrel, limit=this->limit;
    auto l = this->l1;
    int n1 = l[0];
    int u1 = l[1];
    int m1 = l[2];
    l = this->l2;
    int n2 = l[0];
    int u2 = l[1];
    int m2 = l[2];
    l = this->l3;
    int n3 = l[0];
    int u3 = l[1];
    int m3 = l[2];
    l = this->l4;
    int n4 = l[0];
    int u4 = l[1];
    int m4 = l[2];

    if (c1<1e-10||c2<1e-10){ // there are nearly nonexisting coefficients that just bloat the calculations
        return 0;
    }
    //first two quantum numbers follow normalization rules, both longitudinal and transversal wavefunctions don't depend on radial coordinate
    //due to how vectors are generated (see generateCombinations function in variables.cpp) we must only consider the case when eigenfunctions differ in following manner. The same holds in integrate_delta_dispersion_helper_r2
    if (m1!=m3 || m2!=m4){
        return 0;
    }
    if (n1!=n3 || n2!=n4){
        return 0;
    }
    IntegrationWorkspace wsp1(this->limit);
    IntegrationWorkspace wsp2(this->limit);

    // nested integral can be easily turned into 1D integral with bunch of if statements speeding up the code probably 2-fold, though it would be more cumbersome problemshooting-wise. So that I leave for future implementation

    double normsym = 0.5;
    if (m1 == m2 && m2 == m3 && m3 == m4){
        normsym = 0.25;
    }

    auto outer = make_gsl_function( [&](double rho2) {
        auto inner = make_gsl_function( [&](double rho1) {
                        // inside of a function
                        double psi11 = this->radial_function(rho1, n1, u1, m1);
                        double psi12 = this->radial_function(rho2, n1, u1, m1);

                        double psi21 = this->radial_function(rho1, n2, u2, m2);
                        double psi22 = this->radial_function(rho2, n2, u2, m2);

                        double psi31 = this->radial_function(rho1, n3, u3, m3);
                        double psi32 = this->radial_function(rho2, n3, u3, m3);

                        double psi41 = this->radial_function(rho1, n4, u4, m4);
                        double psi42 = this->radial_function(rho2, n4, u4, m4);

                        double val1,val2,val3,val4;

                        if(centerOfMass){
                            val1 = c1*c2*psi11*psi22*psi31*psi42*(rho1+rho2)*0.5;
                            val2 = c1*c2*psi12*psi21*psi32*psi41*(rho1+rho2)*0.5;
                            val3 = c1*c2*psi11*psi22*psi32*psi41*(rho1+rho2)*0.5;
                            val4 = c1*c2*psi12*psi21*psi31*psi42*(rho1+rho2)*0.5;
                        } else{                        
                            val1 = c1*c2*psi11*psi22*psi31*psi42*rho1;
                            val2 = c1*c2*psi12*psi21*psi32*psi41*rho1;
                            val3 = c1*c2*psi11*psi22*psi32*psi41*rho1;
                            val4 = c1*c2*psi12*psi21*psi31*psi42*rho1;
                        }

                        double finval = 0;
                        if(m1==m2 && m3==m4){
                            finval += val1;
                            finval += val2;
                        }
                        if(m1==m3 && m2==m4){
                            finval += val3;
                            finval += val4;
                        }
                        
                        return finval;
                        //return c1*c2*(psi11*psi22+psi12*psi21)*(psi31*psi42+psi31*psi42)*(rho1+rho2)*0.5;
                    } );
        //gsl_integration_qags(inner, 0, 20, epsabs, epsrel, limit, wsp1, &result_inner, &abserr_inner);
                    
        gsl_integration_qag(inner, 0, 20, epsabs, epsrel, limit, this->key, wsp1, &result_inner, &abserr_inner);
                    
        //gsl_integration_qagiu(inner, 0, epsabs, epsrel, limit, wsp1, &result_inner, &abserr_inner);
                                        
        return result_inner;
    });
    //gsl_integration_qags(outer, 0, 20, epsabs, epsrel, limit, wsp2, &result, &abserr);
    
    gsl_integration_qag(outer, 0, 20, epsabs, epsrel, limit, this->key, wsp2, &result, &abserr);
    
    //gsl_integration_qagiu(outer, 0, epsabs, epsrel, limit, wsp1, &result, &abserr);
    

    return result*normsym;
}

double CompleteIntegral::integrate_delta_dispersion_helper_r2(double c1, double c2, bool centerOfMass){
    double result, abserr, result_inner, abserr_inner;
    double epsabs=this->epsabs, epsrel=this->epsrel, limit=this->limit;
    auto l = this->l1;
    int n1 = l[0];
    int u1 = l[1];
    int m1 = l[2];
    l = this->l2;
    int n2 = l[0];
    int u2 = l[1];
    int m2 = l[2];
    l = this->l3;
    int n3 = l[0];
    int u3 = l[1];
    int m3 = l[2];
    l = this->l4;
    int n4 = l[0];
    int u4 = l[1];
    int m4 = l[2];
    if (c1<1e-10||c2<1e-10){ // there are nearly nonexisting coefficients that just bloat the calculations
        return 0;
    }
    //for explanation of following "if" statements implementation read integrate_delta_dispersion_helper_r
    if (m1!=m3 || m2!=m4){
        return 0;
    }
    if (n1!=n3 || n2!=n4){
        return 0;
    }
    IntegrationWorkspace wsp1(this->limit);
    IntegrationWorkspace wsp2(this->limit);


    double normsym = 0.5;
    if (m1 == m2 && m2 == m3 && m3 == m4){
        normsym = 0.25;
    }

    auto outer = make_gsl_function( [&](double rho2) {
        auto inner = make_gsl_function( [&](double rho1) {
                        // inside of a function
                        
                        double psi11 = this->radial_function(rho1, n1, u1, m1);
                        double psi12 = this->radial_function(rho2, n1, u1, m1);

                        double psi21 = this->radial_function(rho1, n2, u2, m2);
                        double psi22 = this->radial_function(rho2, n2, u2, m2);

                        double psi31 = this->radial_function(rho1, n3, u3, m3);
                        double psi32 = this->radial_function(rho2, n3, u3, m3);

                        double psi41 = this->radial_function(rho1, n4, u4, m4);
                        double psi42 = this->radial_function(rho2, n4, u4, m4);

                        double val1, val2, val3, val4;

                        if(centerOfMass){
                            val1 = c1*c2*psi11*psi22*psi31*psi42*(rho1*rho1+rho2*rho2+2*rho1*rho2)*0.25;
                            val2 = c1*c2*psi12*psi21*psi32*psi41*(rho1*rho1+rho2*rho2+2*rho1*rho2)*0.25;
                            val3 = c1*c2*psi11*psi22*psi32*psi41*(rho1*rho1+rho2*rho2+2*rho1*rho2)*0.25;
                            val4 = c1*c2*psi12*psi21*psi31*psi42*(rho1*rho1+rho2*rho2+2*rho1*rho2)*0.25;
                        } else{
                            val1 = c1*c2*psi11*psi22*psi31*psi42*(rho1*rho1);
                            val2 = c1*c2*psi12*psi21*psi32*psi41*(rho1*rho1);
                            val3 = c1*c2*psi11*psi22*psi32*psi41*(rho1*rho1);
                            val4 = c1*c2*psi12*psi21*psi31*psi42*(rho1*rho1);
                        }

                        double finval = 0;
                        if(m1==m2 && m3==m4){
                            finval += val1;
                            finval += val2;
                        }
                        if(m1==m3 && m2==m4){
                            finval += val3;
                            finval += val4;
                        }
                        return finval;
                        //return c1*c2*(psi11*psi22+psi12*psi21)*(psi31*psi41+psi32*psi42)*(rho1*rho1+rho2*rho2+2*rho1*rho2)*0.25;
                    } );
        //gsl_integration_qags(inner, 0, 20, epsabs, epsrel, limit, wsp1, &result_inner, &abserr_inner);
                    
        gsl_integration_qag(inner, 0, 20, epsabs, epsrel, limit, this->key, wsp1, &result_inner, &abserr_inner);
                    
        //gsl_integration_qagiu(inner, 0, epsabs, epsrel, limit, wsp1, &result_inner, &abserr_inner);
                                        
        return result_inner;
    });
    //gsl_integration_qags(outer, 0, 20, epsabs, epsrel, limit, wsp2, &result, &abserr);
    
    gsl_integration_qag(outer, 0, 20, epsabs, epsrel, limit, this->key, wsp2, &result, &abserr);
    
    //gsl_integration_qagiu(outer, 0, epsabs, epsrel, limit, wsp1, &result, &abserr);
    
    return result*normsym;
    
}

