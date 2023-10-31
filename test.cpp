#include <boost/math/special_functions/math_fwd.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <functional>
#include <map>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/factorials.hpp>
#include <tuple>

double e = M_E;
double pi = M_PI;

double particleSize = 1e-2; // to avoid singularity at z=0, temporary approximation

double epsabs = 1e-4;
double epsrel = 1e-4;
double length = 1;
double invlength = 1/length;
size_t limit = 100;

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

double integrateFin(int a, double z, double r2){
    double result, abserr;
    
    double lowBound = abs(z)+particleSize;
    double highBound = std::sqrt(z*z+r2*r2);
    

    IntegrationWorkspace wspFinite(limit);

    auto integrandFinite = make_gsl_function( [&](double u) {

        double exp1 = std::exp(-u*u);
        double exp2 = std::exp(2*r2*std::sqrt(u*u-z*z));
        double rest = std::exp(-2*r2*r2) * std::exp(-z*z);
        double frac1 = (u*u-3*z*z)/(u*u*u*u);
        //change to fast power:
        double frac2 = std::pow((-std::sqrt(u*u-z*z)+r2),a)/std::sqrt(u*u-z*z);
        return exp1*exp2*frac1*frac2*rest;

    });

    gsl_integration_qags(integrandFinite, lowBound, highBound, epsabs, epsrel, limit, wspFinite, &result, &abserr);

    return result;
}

double integrateInf(int a, double z, double r2){
    double result, abserr;
    double result2, abserr2;
    
    double lowBound = abs(z)+particleSize;
    double midBound = abs(z)+1.; // used for integration over singularity before integrating over whole space
    IntegrationWorkspace wspInfinite(limit);

    auto integrandInfinite = make_gsl_function( [&](double u) {

        double exp1 = std::exp(-u*u);
        double exp2 = std::exp(-2*r2*std::sqrt(u*u-z*z));
        double rest = std::exp(-2*r2*r2) * std::exp(-z*z);
        double frac1 = (u*u-3*z*z)/(u*u*u*u);
        //change to fast power:
        double frac2 = std::pow((std::sqrt(u*u-z*z)+r2),a)/std::sqrt(u*u-z*z);
        return exp1*exp2*frac1*frac2*rest;

    });

    // two integrals are used, but in reality only the first one is important if we are considering proper boundaries (integrand goes to 0 pretty fast), thus for better performance we could drop the second one
    gsl_integration_qags(integrandInfinite, lowBound, midBound, epsabs, epsrel, limit, wspInfinite, &result, &abserr);

    gsl_integration_qagiu(integrandInfinite, midBound, epsabs, epsrel, limit, wspInfinite, &result2, &abserr2);

    return result+result2;


}

double integrate_single_term(int a, double z, double r2){
    double resultFin = integrateFin(a, z, r2);
    double resultInf = integrateInf(a, z, r2);
    return resultFin+resultInf;
}

//funkcja do usunięcia
double raising_factorial(double x, int i){
    return boost::math::rising_factorial(x,i);
}

//tutaj powinienem parsować tylko to co potrzebne - wcześniej zbuduję pojedyńczy wektor z polynomiali
double integrate_fourth(double r2, double z, std::vector<std::tuple<double, int>> polynomial){
    double result=0;
    //std::cout << "pol length: " << polynomial.size() << std::endl;
    //teraz należy policzyć całki dla każdego z elementów wielomianu
    //przy kompilacji trzeba zwrócić uwagę na wersję C++17 lub wyższą
    for (const auto& [A, a] : polynomial) {
        result += integrate_single_term(a+1, z, r2); // +1 added because of how scalar product in cylindrical coordinates works
    }
    return result;
}
//całkowanie pozostałych powinno być osiągnięte całką 3 rzędu nestowaną w ten sam sposób co pozostałe.

//to idzie do pliku variables.h i .cpp
std::vector<std::tuple<double, int>> multiplyPolynomial(std::vector<std::tuple<double, int>> pol1, std::vector<std::tuple<double, int>> pol2){
    std::vector<std::tuple<double, int>> resultingPolynomial;

    for(int i=0 ; i<pol1.size() ; ++i){
        for(int j=0 ; j<pol2.size() ; ++j){
            resultingPolynomial.push_back(
                    (std::tuple<double, int>) { //accessing tuple contents by index done with get, two object array casted as tuple, there are better ways of doing that with C++17
                        std::get<0>(pol1.at(i))*std::get<0>(pol2.at(j)),
                        std::get<1>(pol1.at(i))+std::get<1>(pol2.at(j))
                    }
                );
        }
    }
    return resultingPolynomial;
}

double fast_power(double a, int b) {
  double result = 1;
  for (int i = 0; i < b; ++i) {
    result *= a;
  }
  return result;
}

double computePolynomial(double x, std::vector<std::tuple<double, int>> pol){
    double result = 0;
    for(int i=0 ; i<pol.size() ; ++i){
        result += std::get<0>(pol.at(i)) + fast_power(x, std::get<1>(pol.at(i)));
    }
    return result;
}

//liczby kwantowe będą parsowane za pomocą klasy Integral
double integrate_full_dipole(int n1, int u1, int m1, int n2, int u2, int m2, int n3, int u3, int m3, int n4, int u4, int m4){
    //poniższe wielomiany powinny być kombinowane w taki którego nie będę się wstydził postawić wewnątrz czwartej całki
    //=============================
    std::vector<std::tuple<double, int>> polynomial1;
    //generowanie wag i eksponensów wielomianu
    for (int xi = 0 ; xi <= u1 ; ++xi){
        double A = boost::math::rising_factorial(-u1,xi)/boost::math::rising_factorial(abs(m1)+1,xi)/boost::math::factorial<double>(xi);
        std::tuple<double, int> element = {A,int(2*xi)};
        polynomial1.push_back(element);
    }

    std::vector<std::tuple<double, int>> polynomial2;
    for (int xi = 0 ; xi <= u2 ; ++xi){
        double A = boost::math::rising_factorial(-u2,xi)/boost::math::rising_factorial(abs(m2)+1,xi)/boost::math::factorial<double>(xi);
        std::tuple<double, int> element = {A,int(2*xi)};
        polynomial2.push_back(element);
    }

    std::vector<std::tuple<double, int>> polynomial3;
    for (int xi = 0 ; xi <= u3 ; ++xi){
        double A = boost::math::rising_factorial(-u3,xi)/boost::math::rising_factorial(abs(m3)+1,xi)/boost::math::factorial<double>(xi);
        std::tuple<double, int> element = {A,int(2*xi)};
        polynomial3.push_back(element);
    }

    std::vector<std::tuple<double, int>> polynomial4;
    for (int xi = 0 ; xi <= u4 ; ++xi){
        double A = boost::math::rising_factorial(-u4,xi)/boost::math::rising_factorial(abs(m4)+1,xi)/boost::math::factorial<double>(xi);
        std::tuple<double, int> element = {A,int(2*xi)};
        polynomial4.push_back(element);
    }
    //=============================
    std::vector<std::tuple<double, int>> pol_gamma = multiplyPolynomial(polynomial1, polynomial2);
    std::vector<std::tuple<double, int>> pol_a = multiplyPolynomial(polynomial3, polynomial4);
    //std::vector<std::tuple<double, int>> polynomial = multiplyPolynomial(pol_a, pol_gamma);

    std::cout << "polynomial gamma length: " << pol_gamma.size() <<"\npolynomial a length: " << pol_a.size() << std::endl;

    double result, abserr, resultInner, abserrInner, resultMiddle, abserrMiddle, resultInner2, abserrInner2;

    // niepewności być może trzeba będzie zmienić żeby zaoszczędzić czas liczenia

    //double particleSize = 0;

    IntegrationWorkspace wspInner(limit);
    IntegrationWorkspace wspMiddle(limit);
    IntegrationWorkspace wspOuter(limit);

    int temp_counter = 0;

    auto outer = make_gsl_function( [&](double r1) {
        auto middle = make_gsl_function( [&](double z1) {
            auto inner = make_gsl_function( [&](double boldZ) {

                //to oczywiście nie jest pełna forma, jak wcześniej muszę pamiętać o symetryzacji i policzyć cztery różne fragmenty całki
                //poza tym będę musiał dodać delty związane z liczbami momentopędowymi
                //ponadto muszę pamiętać o części zespolonej, a przynajmniej sprawdzeniu, czy faktycznie dąży do zera
                // do cosinusa trzeba dodać sinus dla spełnienia założeń zespolonych
                double Zpart = std::cos(2*pi*invlength * (-n1*z1 + n3*z1 - n2*(boldZ+z1) + n4*(boldZ+z1)));
                //nie wiem czy compute polynomial nie wywołuje się jednokrotnie
                //wszystko powinienem poprzekładać w swoim czasie na zewnątrz tych całek dla poprawy wydajności
                //exp(-2*r1*r1) is inside fourth integral
                double r1part = computePolynomial(r1, pol_gamma);
                //double r2Leftovers = std::exp(-(boldZ)*(boldZ));
                double r2Leftovers = 1;
                double r2part = integrate_fourth(r1, boldZ, pol_a);
                ++temp_counter;
                //std::cout << r2part << std::endl;
                //std::cout << "bold Z value : " << boldZ << "\tintegral value : " << Zpart*r1part*r2Leftovers*r2part << std::endl;
                return Zpart*r1part*r2Leftovers*r2part;
            });
            //araising singularities should vanish after dropping region close to 0 or solving analitically fourth integral for specific case (last operation lets us stop gsl algorithm from derailing)
            //std::cout << "z1 value : " << z1 << std::endl;
            gsl_integration_qags(inner, -z1, length-z1, epsabs, epsrel, limit, wspInner, &resultInner, &abserrInner);
            //std::cout << "value of z1: " << z1 << " result: " << resultInner+resultInner2 << std::endl;
            return resultInner;
        });
        std::cout << "r2 value:\t" << r1 << std::endl;
        gsl_integration_qags(middle, 0, length, epsabs, epsrel, limit, wspMiddle, &resultMiddle, &abserrMiddle);
        std::cout << "\tresult of 3'rd rank integral with r1=" << r1 << " : " << resultMiddle << std::endl;
        return resultMiddle;
    });
    std::cout << "first" << std::endl;
    //gsl_integration_qagiu(outer,0 , epsabs, epsrel, limit, wspOuter, &result, &abserr);
    gsl_integration_qags(outer, particleSize, 10 , epsabs, epsrel, limit, wspOuter, &result, &abserr);
    std::cout << "second" << std::endl;
    std::cout << "results: " << resultInner << ' ' << resultMiddle << ' ' << result << std::endl;
    std::cout << "counter: " << temp_counter << std::endl;
    return result;
}

int main(){

    int a = 1;
    double z = 1;
    double r2 = 2;

    double integralfin = integrateFin(a, z, r2);
    double integralinf = integrateInf(a, z, r2);
    
    std::cout << "Finite integral result: " << integralfin << "\nInfinite integral reult: " << integralinf << std::endl;

    int n1=1;
    int u1=1;
    int m1=1;
    int n2=1;
    int u2=1;
    int m2=1;
    int n3=1;
    int u3=1;
    int m3=1;
    int n4=1;
    int u4=1;
    int m4=1;

    double integral_full = integrate_full_dipole(n1,u1,m1,n2,u2,m2,n3,u3,m3,n4,u4,m4);

    std::cout << "result of fourth order integral: " << integral_full << std::endl;

    return 0;
}