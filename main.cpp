#include "functions.h"
#include <stdio.h>  // printf, sprintf, fopen, fprintf, fclose
#include <iostream> // cout, endl, ...
#include <cmath>    // pow
#include <ctime>    // clock, CLOCKS_PER_SEC
#include <string>   // string, stod

using namespace std;

void load_input(int &R, int &n_start, int &n_gen, int &D, int &NU_X, int &NU_alpha, int &n_simulations, vector<double> &eps, vector<double> &s, int argc, char *argv[]);

int main(int argc, char *argv[]){

    if (argc < 10){cerr << "ERROR: the arguments aren't enough !\n"; return -1;}

    // --------------- parameters --------------- //
    // we have R resources and an ancestral population consisting
    // of n_start individuals and we study them for n_gen generations
    int R; // number of resources
    int n_start; // initial number of individuals
    int n_gen; // number of generations
    int D; // depth of sequencing

    // number of mutations in the population per generation
    int NU_X;
    int NU_alpha;

    int n_simulations; // number of simulations with the same parameters

    vector<double> eps; // small perturbations around the completely symmetric state for beta
    vector<double> s; // fitness increment in mutation

    // we load the input from command line
    load_input(R, n_start, n_gen, D, NU_X, NU_alpha, n_simulations, eps, s, argc, argv);
    
    output_run out;
    
    clock_t tstart; // used to time the performances
    char filename[50];
    FILE * fp;

    // --------------- run --------------- //
    for (size_t i = 0; i < eps.size(); i++){

        sprintf(filename, "./Data/n_species_%.1e_%.1e.txt", eps[i], s[i]);
        fp = fopen(filename, "w");
        
        tstart = clock();
        cout << "----------------------------------------\n";
        printf("eps = %.1e, s = %.1e\t(%li/%li)\n", eps[i], s[i], i + 1, eps.size());

        for (int j = 0; j < n_simulations; j++){

            printf("  simulation #%i out of %i\n", j + 1, n_simulations);
            printf("  time since the beginning for these values = %f s\n", (clock() - tstart) / (double) CLOCKS_PER_SEC);

            if(j > 0){
                printf("  Results for the previous simulation:\n");
                printf("    # of sequenced species = %i\n", out.n_species);
            }

            out = run(R, n_start, n_gen, D, NU_X, NU_alpha, eps[i], s[i]);

            fprintf(fp, "%f, %i\n", out.scaled_rate, out.n_species);

            cout << endl;

        }

        fclose(fp);

    }

    return 0;

}

// function to read arguments from command line
void load_input(int &R, int &n_start, int &n_gen, int &D, int &NU_X, int &NU_alpha, int &n_simulations, vector<double> &eps, vector<double> &s, int argc, char *argv[]){

    R = int(stod(argv[1])); // number of resources
    n_start = int(stod(argv[2])); // initial number of individuals
    n_gen = int(stod(argv[3])); // number of generations
    D = int(stod(argv[4])); // depth of sequencing

    // number of mutations in the population per generation
    NU_X = int(stod(argv[5]));
    NU_alpha = int(stod(argv[6]));

    n_simulations = int(stod(argv[7])); // number of simulations with the same parameters

    vector<string> data(argv + 8, argv + argc + !argc); // load the remaining data in this vector

    // split data between eps and s
    std::transform(data.begin(), data.end() - data.size() / 2, std::back_inserter(eps),
               [](const std::string& str) { return std::stod(str); });

    std::transform(data.begin() + data.size() / 2, data.end(), std::back_inserter(s),
               [](const std::string& str) { return std::stod(str); });

}