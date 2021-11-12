// ------------------------------------------------------------- //
// *** WARNING ***: to compile this file make sure to change the //
// LIBS value in the Makefile to the path of your BOOST library! //
// ------------------------------------------------------------- //

#include "functions.h"
#include <stdio.h>  // printf, sprintf, fopen, fprintf, fclose
#include <iostream> // cout, endl, ...
#include <cmath>    // pow
#include <ctime>    // clock, CLOCKS_PER_SEC
#include <string>   // string, stod

using namespace std;

void load_input(int &R, int &n_start, int &n_gen, int &dt_gen, int &D, double &NU_Y, double &NU_alpha, int &n_simulations, int &fast, int &generalist, int &random_selection, int &seed, double &chi, double &eps, double &s, int argc, char *argv[]);

int main(int argc, char *argv[]){

    if (argc < 15){cerr << "\nERROR: the arguments aren't enough !\n\n"; return -1;}

    // --------------- parameters --------------- //
    // we have R resources and an ancestral population consisting
    // of n_start individuals and we study them for n_gen generations
    int R; // number of resources
    int n_start; // initial number of individuals
    int n_gen; // number of generations

    int dt_gen; // number of generations between two sequencings
    int D; // depth of sequencing

    // number of mutations in the population per generation
    double NU_Y; // number of fitness mutations
    double NU_alpha; // number of strategy mutations

    int n_simulations; // number of simulations with the same parameters

    int fast; // if we want faster simulations we have to give up "checkpoints" (sequencing every dt_gen generations) (fast = 1, slow = 0)
    int generalist; // specialist = 0, generalist = 1, complete pool = 2
    int random_selection; // should the selection step be stochastic or deterministic? (stochastic = 1, deterministic = 0)
    int seed; // seed for the random number generator

    double chi; // cost per resource (death rate)

    double eps; // small perturbations around beta
    double s; // fitness increment in mutation

    // we load the input from command line
    load_input(R, n_start, n_gen, dt_gen, D, NU_Y, NU_alpha, n_simulations, fast, generalist, random_selection, seed, chi, eps, s, argc, argv);
    
    output_run out;
    
    clock_t tstart; // used to time the performances
    char filename[100];
    FILE * fp;

    // --------------- run --------------- //

    // initialize the seed of the random number generator
    RNG::rng.seed(seed);

    tstart = clock();
    cout << "------------------------------------\n";
    printf("eps = %.1e, s = %.1e\n", eps, s);
    printf("NU_Y = %.1e, NU_alpha = %.1e\n", NU_Y, NU_alpha);
    printf("chi = %.2f\n\n", chi);

    for (int j = 0; j < n_simulations; j++){

        printf("  simulation #%i out of %i\n", j + 1, n_simulations);
        
        out = run(R, n_start, n_gen, dt_gen, D, NU_Y, NU_alpha, eps, s, fast, generalist, random_selection, chi); // the simulations happen here!
        
        sprintf(filename, "./Data/run_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e.txt", (double) R, (double) n_start, (double) n_gen, (double) generalist, (double) random_selection, eps, s, NU_Y, NU_alpha, chi);
        fp = fopen(filename, "a");

        fprintf(fp, "%f\n", out.scaled_rate);
        fprintf(fp, "%i\n", out.n_species);

        fclose(fp);

        printf("  time for these values = %.2f s\n", (clock() - tstart) / (double) CLOCKS_PER_SEC);
        cout << endl;

    }

    printf("Total time required = %.2f s\n", (clock() - tstart) / (double) CLOCKS_PER_SEC);

    return 0;

}

// function to read arguments from command line
void load_input(int &R, int &n_start, int &n_gen, int &dt_gen, int &D, double &NU_Y, double &NU_alpha, int &n_simulations, int &fast, int &generalist, int &random_selection, int &seed, double &chi, double &eps, double &s, int argc, char *argv[]){

    R = int(stod(argv[1])); // number of resources
    n_start = int(stod(argv[2])); // initial number of individuals
    n_gen = int(stod(argv[3])); // number of generations

    dt_gen = int(stod(argv[4])); // number of generations between two sequencings
    D = int(stod(argv[5])); // depth of sequencing

    // number of mutations in the population per generation
    NU_Y = stod(argv[6]);
    NU_alpha = stod(argv[7]);

    n_simulations = int(stod(argv[8])); // number of simulations with the same parameters

    fast = int(stod(argv[9])); // if we want faster simulations we have to give up "checkpoints" (sequencing every dt_gen generations) (fast = 1, slow = 0)    
    generalist = int(stod(argv[10])); // specialist = 0, generalist = 1, complete pool = 2
    random_selection = int(stod(argv[11])); // should the selection step be stochastic or deterministic? (stochastic = 1, deterministic = 0)
    seed = int(stod(argv[12])); // seed for the random number generator

    chi = stod(argv[13]); // cost per resource (death rate)

    eps = stod(argv[14]); // small perturbations around beta
    s = stod(argv[15]); // fitness increment in mutation

}