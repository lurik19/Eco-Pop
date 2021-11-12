// ------------------------------------------------------------- //
// *** WARNING ***: to compile this file make sure to change the //
// LIBS value in the Makefile to the path of your BOOST library! //
// ------------------------------------------------------------- //

#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <vector>                       // vector
#include <map>                          // map
#include <numeric>                      // accumulate, partial_sum, reduce
#include <algorithm>                    // transform, copy, count, max
#include <cmath>                        // exp, log
#include <stdio.h>                      // sprintf, fopen, fprintf, fclose
#include <iostream>                     // cout, ...

// BOOST LIBRARIES
#include <boost/functional/hash.hpp>            // boost::hash_range
#include <boost/random/mersenne_twister.hpp>    // boost::random::mt19937_64
#include <boost/unordered_map.hpp>              // boost::unordered_map
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>


using namespace std;

struct output_run{
    double scaled_rate;
    int n_species;
};

// hashing implemented to use vectors as keys in the unordered_map
template <typename Container>
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};

// data structure that we will use a lot in the program
typedef boost::unordered_map<vector<double>, double, container_hash<vector<double>>> u_map;

namespace RNG{
    extern boost::random::mt19937_64 rng; // not sure if this is orthodox C++
}

// -------------------- Definitions of functions --------------------//
// functions that initialize the system
void initial_alpha_Y_gen(u_map &u_map_alpha_Y, int R, int n);
void initial_alpha_Y_spec(u_map &u_map_alpha_Y, int R, int n, double s);
void complete_pool(u_map &u_map_alpha_Y, int R, int n_start, double chi);
int binomialCoefficients(int n, int k);
vector<double> compute_beta(int R, double eps);

// ecological selection
void selection_step(u_map &u_map_alpha_Y, vector<double> &beta, int R, int n_start, int random_selection, double chi);

// evolution: mutations
void mutation_step(u_map &u_map_alpha_Y, double s, double U_Y, double U_alpha, int R);
vector<vector<double>> fitness_mutation(vector<double> alpha_Y_mu, int k_Y, double s, int R);
vector<vector<double>> lose_mutation(vector<double> alpha_Y_mu, int k_lose, int R);
vector<vector<double>> gain_mutation(vector<double> alpha_Y_mu, int k_gain, int R);
void update_umap(const vector<vector<double>> &alpha_Y_new, u_map &u_map_alpha_Y);

// extract infos from unordered_map
vector<vector<double>> get_alpha_Y(u_map &u_map_alpha_Y);
vector<vector<double>> get_alpha(u_map &u_map_alpha_Y, int R);
vector<double> get_Y(u_map &u_map_alpha_Y, int R);
vector<int> get_n(u_map &u_map_alpha_Y);

// compute Lyapunov function
double lyapunov_function(u_map &u_map_alpha_Y, double chi, int R, vector<double> beta);

// summary: write on files
int summary(u_map u_map_alpha_Y, int D, int R, int n_start, int n_gen, int generalist, int random_selection, int write_on_file, double NU_Y, double NU_alpha, double eps, double s, double chi, vector<double> beta);

// run the system
output_run run(int R, int n_start, int n_gen, int dt_gen, int D, double NU_Y, double NU_alpha, double eps, double s, int fast, int generalist, int random_selection, double chi);
vector<double> compute_eY_i(u_map &u_map_alpha_Y, vector<double> &beta, int R);
double std_dev(vector<double> &v);

// print data structures
void print_umap(u_map &umap);
void print_map(map<int, vector<double>> map);
template <typename T>
void print_vec(vector<T> v);

# endif