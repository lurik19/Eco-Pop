// ------------------------------------------------------------- //
// *** WARNING ***: to compile this file make sure to change the //
// LIBS value in the Makefile to the path of your BOOST library! //
// ------------------------------------------------------------- //

#ifndef __FUNCTIONS__
#define __FUNCTIONS__
#include <random>       // default_random_engine, normal_distribution, poisson_distribution, binomial_distribution, uniform_int_distribution
#include <vector>       // vector
#include <numeric>      // accumulate, partial_sum
#include <algorithm>    // transform, copy, count
#include <map>          // map
#include <unordered_map>          // unordered_map
#include <cmath>        // exp, log
#include <boost/functional/hash.hpp> // hash_range

#include <iostream>
#include <stdio.h>  // sprintf, fopen, fprintf, fclose


using namespace std;

struct output_run{
    double scaled_rate;
    int n_species;
};

// hashing to use vectors as keys in the unordered_map
template <typename Container>
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};

typedef unordered_map<vector<double>, int, container_hash<vector<double>>> u_map;

// -------------------- Definitions of functions --------------------//
vector<double> initial_alpha_X(int R);
vector<double> compute_beta(int R, double eps);

void selection_step(u_map &u_map_alpha_X, vector<double> &beta, int R);
vector<vector<double>> get_alpha_X(u_map &u_map_alpha_X, int R);
vector<vector<double>> get_alpha(u_map &u_map_alpha_X, int R);
vector<double> get_X(u_map &u_map_alpha_X, int R);
vector<int> get_n(u_map &u_map_alpha_X);

void mutation_step(u_map &u_map_alpha_X, double s, double U_X, double U_alpha, int R);
vector<vector<double>> fitness_mutation(vector<double> alpha_X_mu, int k_X, double s, int R);
vector<vector<double>> lose_mutation(vector<double> alpha_X_mu, int k_lose, int R);
vector<vector<double>> gain_mutation(vector<double> &alpha_X_mu, int k_gain, int R);
void update_umap(vector<vector<double>> &alpha_X_new, u_map &u_map_alpha_X);

int sequencing(u_map u_map_alpha_X, int D, int R);

output_run run(int R, int n, int n_gen, int D, int NU_X, int NU_alpha, double eps, double s);
output_run run_complete(int R, int n, int n_gen, int dt_gen, int D, int NU_X, int NU_alpha, double eps, double s);

# endif