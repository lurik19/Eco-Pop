// ------------------------------------------------------------- //
// *** WARNING ***: to compile this file make sure to change the //
// LIBS value in the Makefile to the path of your BOOST library! //
// ------------------------------------------------------------- //

#include "functions.h"

default_random_engine rng(0); // not sure if this is orthodox C++

// -------------------- Implementations of functions --------------------//

vector<double> initial_alpha_X(int R){

    vector<double> alpha_X (R + 1, 1. / R); //generalist strategy
    int X = 0; // the fitness value of the first strain does not influence the results
               // (X_1 always cancels in \lambda_\mu, because \overline{X}_i = X_1 + log(...))

    alpha_X[R] = X;

    return alpha_X;

}

vector<double> compute_beta(int R, double eps){

    vector<double> epss (R);
    vector<double> beta (R);
    normal_distribution<double> dist(0, eps);

    for (int i = 0; i < R; i++){

        epss [i] = dist(rng);

    }
    
    double sum = accumulate(epss.begin(), epss.end(), decltype(epss)::value_type(0));
    transform(epss.begin(), epss.end(), beta.begin(), [R, sum](double &v){ return 1. / R * (1 + v - 1. / R * sum);}); // Eq (S92)

    return beta;

}

void selection_step(u_map &u_map_alpha_X, vector<double> &beta, int R){

    // extract infos from u_map_alpha_X
    vector<vector<double>> alpha_X = get_alpha_X(u_map_alpha_X, R);
    vector<int> n = get_n(u_map_alpha_X);
    
    int n_size = n.size();

    vector<double> X_i(R, 0);
    vector<double> lambda(n_size, 0);

    int sum_n = accumulate(n.begin(), n.end(), decltype(n)::value_type(0));

    // let's compute \overline{X}_i(t)
    for (int i = 0; i < R; i++){
        for (int mu = 0; mu < n_size; mu++){
            X_i[i] += alpha_X[mu][i] * exp(alpha_X[mu][R]) * n[mu]; // X_i[i] += alpha[mu][i] * exp(X[mu]) * n[mu]
        }
        X_i[i] = log(X_i[i] / (beta[i] * sum_n));
    }

    // let's compute \lambda_\mu
    for (int mu = 0; mu < n_size; mu++){
        for (int i = 0; i < R; i++){
            lambda[mu] += alpha_X[mu][i] * exp(alpha_X[mu][R] - X_i[i]) * n[mu]; // lambda[mu] += alpha[mu][i] * exp(X[mu] - X_i[i]) * n[mu]
        }

        // let's sample new strain sizes from a poisson distribution
        poisson_distribution<int> dist(lambda[mu]);
        n[mu] = dist(rng);
    }

    int mu = 0;
    for(auto& imap: u_map_alpha_X){

        imap.second = n[mu];
        mu ++;

    }

}

vector<vector<double>> get_alpha_X(u_map &u_map_alpha_X, int R){

    vector<vector<double>> alpha_X(u_map_alpha_X.size(), vector<double>(R));

    size_t i = 0;
    for(auto const& imap: u_map_alpha_X){
        
        alpha_X[i] = imap.first;
        i++;
    
    }

    return alpha_X;

}

vector<vector<double>> get_alpha(u_map &u_map_alpha_X, int R){

    vector<vector<double>> alpha_X = get_alpha_X(u_map_alpha_X, R);
    vector<vector<double>> alpha(alpha_X.size(), vector<double>(R));
    
    for (size_t mu = 0; mu < alpha_X.size(); mu++)
        copy(alpha_X[mu].begin(), alpha_X[mu].end() - 1, alpha[mu].begin());

    return alpha;

}

vector<double> get_X(u_map &u_map_alpha_X, int R){

    vector<vector<double>> alpha_X = get_alpha_X(u_map_alpha_X, R);
    vector<double> X(alpha_X.size());
    int last = alpha_X[0].size() - 1;

    for (size_t mu = 0; mu < alpha_X.size(); mu++)
        X[mu] = alpha_X[mu][last];

    return X;

}

vector<int> get_n(u_map &u_map_alpha_X){

    vector<int> n(u_map_alpha_X.size());
    int i = 0;
    for(auto const& imap: u_map_alpha_X){
        n[i] = imap.second;
        i++;
    }
    return n;

}

void mutation_step(u_map &u_map_alpha_X, double s, double U_X, double U_alpha, int R){

    // extract infos from map_alpha_X
    vector<vector<double>> alpha_X = get_alpha_X(u_map_alpha_X, R);
    vector<int> n = get_n(u_map_alpha_X);

    int k_X, k_lose, k_gain;

    for (size_t mu = 0; mu < alpha_X.size(); mu++){

        // fitness mutations
        binomial_distribution<int> dist1(n[mu], U_X);
        k_X = dist1(rng);

        // if k_X = 0 it doesn't add anything to u_map_alpha_X (same for strategy mutations)
        vector<vector<double>> alpha_X_new1 = fitness_mutation(alpha_X[mu], k_X, s, R);
        update_umap(alpha_X_new1, u_map_alpha_X);

        // strategy mutations
        // lose
        binomial_distribution<int> dist2(n[mu], U_alpha * (R - count(alpha_X[mu].begin(), alpha_X[mu].end() - 1, 0)));
        k_lose = dist2(rng);

        vector<vector<double>> alpha_X_new2 = lose_mutation(alpha_X[mu], k_lose, R);
        update_umap(alpha_X_new2, u_map_alpha_X);

        // gain
        binomial_distribution<int> dist3(n[mu], U_alpha * R);
        k_gain = dist3(rng);

        vector<vector<double>> alpha_X_new3 = gain_mutation(alpha_X[mu], k_gain, R);
        update_umap(alpha_X_new3, u_map_alpha_X);

        u_map_alpha_X[alpha_X[mu]] -= (k_X + k_gain + k_lose); // we subtract the number of individuals with mutations from the original strain

        if(u_map_alpha_X[alpha_X[mu]] == 0)   // if a strain is extinct
            u_map_alpha_X.erase(alpha_X[mu]); // we remove it
    }

}

vector<vector<double>>fitness_mutation(vector<double> alpha_X_mu, int k_X, double s, int R){

    vector<vector<double>> alpha_X_new(k_X, vector<double>(R + 1));

    alpha_X_mu[R] += s; // we increase the (log)fitness by s

    for (int mu = 0; mu < k_X; mu++){
        alpha_X_new[mu] = alpha_X_mu;
    }

    return alpha_X_new;

}

vector<vector<double>> lose_mutation(vector<double> alpha_X_mu, int k_lose, int R){

    vector<vector<double>> alpha_X_new(k_lose, vector<double>(R + 1));
    vector<int> non_zeros(R + 1);

    int nz_size = 0;
    // we find the positions of the non-zero elements ** in alpha **
    for(int i = 0; i < R; i++){
        if(alpha_X_mu[i]){ // alpha_X_mu[i] != 0
            non_zeros[nz_size] = i;
            nz_size++;
        }
    }

    if(nz_size > 1){ // a strain must consume at least one resource

        uniform_int_distribution<int> dist(0, nz_size - 1);
        vector<double> alpha_X_copy(R + 1);

        // let's normalize the resources consumption
        for(int i = 0; i < R; i++){
            alpha_X_mu[i] *= (double) nz_size / (nz_size - 1);
        }

        for (int mu = 0; mu < k_lose; mu++){
            
            alpha_X_copy = alpha_X_mu;
            // select one resource with I_{\mu, i} = 1 and put I_{\mu, i} = 0
            alpha_X_copy[non_zeros[dist(rng)]] = 0;
            alpha_X_new[mu] = alpha_X_copy;

        }
    }
    
    return alpha_X_new;

}

vector<vector<double>> gain_mutation(vector<double> &alpha_X_mu, int k_gain, int R){

    vector<vector<double>> alpha_X_new(k_gain, vector<double>(R + 1));
    int num_nz_old, num_nz;

    num_nz_old = R - count(alpha_X_mu.begin(), alpha_X_mu.end() - 1, 0); // number of non-zero elements in ** alpha[mu] **

    if (num_nz_old < R){ // this slightly improves the performances

        uniform_int_distribution<int> dist(0, R - 1);
        vector<double> alpha_X_copy(R + 1);

        for (int mu = 0; mu < k_gain; mu++){

            alpha_X_copy = alpha_X_mu;

            // select one resource and put I_{\mu, i} = 1
            alpha_X_copy[dist(rng)] = 1. / num_nz_old;

            num_nz = R - count(alpha_X_copy.begin(), alpha_X_copy.end() - 1, 0); // number of non-zero elements in the new ** alpha **

            // let's normalize the resources consumption
            for(int i = 0; i < R; i++){
                alpha_X_copy[i] *= (double) num_nz_old / num_nz;
            }

            alpha_X_new[mu] = alpha_X_copy;

        }
    } else { // num_nz_old = R

        for (int mu = 0; mu < k_gain; mu++)
            alpha_X_new[mu] = alpha_X_mu;

    }

    return alpha_X_new;

}

void update_umap(vector<vector<double>> &alpha_X_new, u_map &u_map_alpha_X){
    
    for(size_t mu = 0; mu < alpha_X_new.size(); mu++)
        u_map_alpha_X[alpha_X_new[mu]] ++;

}

int sequencing(u_map u_map_alpha_X, int D, int R){

    u_map u_map_count_species;
    vector<vector<double>> alpha = get_alpha(u_map_alpha_X, R);
    vector<int> n = get_n(u_map_alpha_X);

    // let's create the "inverse" map, with the cumulative sum of n's as keys
    partial_sum(n.begin(), n.end(), n.begin());
    map<int, vector<double>> map_n;

    for (size_t i = 0; i < n.size(); i++)
        map_n[n[i]] = alpha[i];

    // let's generate D random numbers between 1 and \sum n_\mu
    vector<int> rands(n[n.size() - 1], 1);
    partial_sum(rands.begin(), rands.end(), rands.begin());
    shuffle(rands.begin(), rands.end(), rng);

    map<int, vector<double>>::iterator idx;

    // let's list the occurences of the sampled ecotypes
    for(int i = 0; i < D; i++){
        idx = map_n.upper_bound(rands[i] - 1);
        u_map_count_species[idx->second] ++;
    }

    // we won't include mutations with only one individual
    auto imap = u_map_count_species.begin();
    while (imap != u_map_count_species.end()){
        if (imap->second == 1){
            // erase() will return the next iterator. We don't need to increment
            imap = u_map_count_species.erase(imap);
        } else{
            // go to the next entry
            imap++;
        }
    }

    return u_map_count_species.size();
    
}

output_run run(int R, int n, int n_gen, int D, int NU_X, int NU_alpha, double eps, double s){

    vector<double> alpha_X = initial_alpha_X(R); // initial resource strategy for the ancestral population
    vector<double> beta = compute_beta(R, eps); // resources' supply rates
    
    u_map u_map_alpha_X {{alpha_X, n}};

    // probabilities of mutations
    double U_X = (double) NU_X / n;
    double U_alpha = (double) NU_alpha / (n * R);

    for (int i = 0; i < n_gen; i++){
        selection_step(u_map_alpha_X, beta, R);
        mutation_step(u_map_alpha_X, s, U_X, U_alpha, R); 
    }

    int n_sequenced_species = sequencing(u_map_alpha_X, D, R);

    output_run result;
    result.n_species = n_sequenced_species;
    result.scaled_rate = U_X * s * s * (double)R / (U_alpha * eps * eps);
    
    return result;

}