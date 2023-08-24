#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <array>

class CompleteIntegral{
    private:
        //int l1[3], l2[3], l3[3], l4[3];
        std::array<int,3> l1,l2,l3,l4; //index triplets - l1 and l2 belong to first vector, l3 and l4 to second vector
        int limit;
        double result, error, epsabs, epsrel;
        double length, rev_length;
        double norm1, norm2, norm3, norm4;// sqrt(1/integral)
        int key;
        static void check_if_valid_key(int k);
        static double fast_power(double,int);
        double dipole_integral_function(std::array<int,3>, std::array<int,3>, std::array<int,3>, std::array<int,3>);
    public:
        CompleteIntegral(int limit, double epsabs, double epsrel){
            this->epsabs = epsabs;
            this->epsrel = epsrel;
            this->limit = limit;
            this->l1 = {0,0,0};
            this->l2 = {0,0,0};
            this->l3 = {0,0,0};
            this->l4 = {0,0,0};
            this->result = 0;
            this->error = 0;
            this->norm1 = 0;
            this->norm2 = 0;
            this->length = 1;
            this->rev_length = 1/this->length;
            this->key = 1;
            this->check_if_valid_key(this->key);
        }
        double get_rev_length(){
            return this->rev_length;
        }
        void change_key(int k){
            this->key = k;
            this->check_if_valid_key(this->key);
        }
        void add_to_key(){
            this->key += 1;
            this->check_if_valid_key(this->key);
        }
        void change_length(double l){
            this->length = l;
            this-> rev_length = 1/l;
        }
        ~CompleteIntegral(){}
        double get_result(){return this->result;}
        void change_norm(double norm1, double norm2, double norm3, double norm4){
            this->norm1 = norm1;
            this->norm2 = norm2;
            this->norm3 = norm3;
            this->norm4 = norm4;
        };
        double integrate_norm_external(int, int, int); // gives integration value for scalar product of single base function
        void changeMultiIndex(int, int, int, int);
        void changeLimit(int);
        void restart_counter(){
            this->result = 0;
            this->error = 0;
        }
        void print_result();

        double fastLaplacianWithHarmonicEnergy(unsigned int, unsigned int);//call only for diagonal elements! arguments are numbers of basis functions
        double integrate_over_delta(double);
        double integrate_over_dipole(double);
        

        double fast_add_over_harmonic(unsigned int, unsigned int);
        double test_integral();
};


#endif /* INTEGRAL_H */