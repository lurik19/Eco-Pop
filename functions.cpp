// ------------------------------------------------------------- //
// *** WARNING ***: to compile this file make sure to change the //
// LIBS value in the Makefile to the path of your BOOST library! //
// ------------------------------------------------------------- //

#include "functions.h" 

// let's define a random number generator (we initialize the seed in main())
// (mt = mersenne twister: this is a reliable algorithm with GOOD RANDOMNESS)
namespace RNG{
    boost::random::mt19937_64 rng; // not sure if this is orthodox C++
}

// -------------------- Implementations of functions -------------------- //

// -------------------- functions that initialize the system -------------------- //
void initial_alpha_Y_gen(u_map &u_map_alpha_Y, int R, int n){
    
    vector<double> alpha_Y (R + 1, 1.); //generalist strategy

    // the fitness value of the first strain does not influence the results
    // (Y_1 always cancels in \lambda_\mu, because \overline{Y}_i = Y_1 + log(...))
    alpha_Y[R] = 0;

    u_map_alpha_Y[alpha_Y] = n;

}

void initial_alpha_Y_spec(u_map &u_map_alpha_Y, int R, int n, double s){

    vector<vector<double>> alpha_Y(R, vector<double>(R + 1, 0));
    boost::random::normal_distribution<double> dist(0, s); // distribution from which we extract fitnesses of the strains

    for(int i = 0; i < R; i++){

        alpha_Y[i][i] = 1; // \alpha_{\mu,i} = \delta_{\mu,i}
        alpha_Y[i][R] = dist(RNG::rng); // fitness

    }

    for (auto const &aY: alpha_Y) 
        u_map_alpha_Y[aY] = n / R;

}

void complete_pool(u_map &u_map_alpha_Y, int R, int n_start, double chi){

    vector<double> zeros(R, 0);
    vector<double> v;
    vector<double> vv;

    for(int i = 1; i <= R; i++){
        v = zeros;
        for(int j = 0; j < i; j++){
            v[j] = 1;
        }
        
        for(int j = 0; j < binomialCoefficients(R, i); j++){
            next_permutation(v.begin(), v.end());

            vv = v;
            vv.push_back(0); // fitness
            u_map_alpha_Y[vv] = n_start / (pow(2, R) - 1);
        }
        
    }

}

int binomialCoefficients(int n, int k) { // this function computes n \choose k
    
    if (k == 0 || k == n)
        return 1;

    return binomialCoefficients(n - 1, k - 1) + binomialCoefficients(n - 1, k);

}

vector<double> compute_beta(int R, double eps){

    vector<double> beta(R, 1);
    partial_sum(beta.begin(), beta.end(), beta.begin());

    // normalization for beta
    double norm = accumulate(beta.begin(), beta.end(), 0.0);
    for(int i = 0; i < R; i++)
        beta[i] /= norm;

    // noise on beta
    vector<double> epss (R);
    boost::random::normal_distribution<double> dist(0, eps);
    for (int i = 0; i < R; i++){
        epss[i] = dist(RNG::rng);
    }
    
    sort(epss.begin(), epss.end()); // we sort the vector in ascending order

    // normalization for the noise on beta
    double sum = 0;
    for(int i = 0; i < R; i++){
        sum += epss[i] * beta[i];
    }

    // let's apply the noise on beta
    for (int i = 0; i < R; i++){
        beta[i] = beta[i] * (1 + epss[i] - sum);
    }

    return beta;

}

// -------------------- ecological selection -------------------- //
void selection_step(u_map &u_map_alpha_Y, vector<double> &beta, int R, int n_start, int random_selection, double chi){

    vector<double> lambda(u_map_alpha_Y.size(), 0);
   
    double eY_i_mu;
    double small = 1E-11; // small number used to avoid divergence

    int mu = 0;
    // let's compute \lambda_\mu
    // lambda[mu] = \sum_i alpha[mu][i] * exp(Y[mu] - Y_i) * n[mu]
    for(auto const& imap_mu: u_map_alpha_Y){
        for (int i = 0; i < R; i++){
            eY_i_mu = 0;
            // eY_i_mu = exp(Y_i - Y[mu]) = \sum_nu alpha[nu][i] * exp(Y[nu] - Y[mu]) * n[nu]
            for(auto const& imap_nu: u_map_alpha_Y)
                eY_i_mu += imap_nu.first[i] * exp(imap_nu.first[R] - imap_mu.first[R]) * imap_nu.second;
            eY_i_mu /= beta[i];
            if(eY_i_mu == 0) eY_i_mu = small;

            lambda[mu] += imap_mu.first[i] / eY_i_mu;
        }
        
        lambda[mu] *= imap_mu.second;
        mu++;

    }

    // normalization which ensures that the total population size remains near n_start ± O(sqrt(n_start))
    double C = n_start / accumulate(lambda.begin(), lambda.end(), decltype(lambda)::value_type(0));

    mu = 0;
    for(auto const& imap: u_map_alpha_Y){
        lambda[mu] *= C;
        // cost per resource (death rate)
        lambda[mu] /= (1 + chi * accumulate(imap.first.begin(), imap.first.end() - 1, 0.));
        mu++;
    }

    mu = 0;
    auto imap = u_map_alpha_Y.begin();

    while (imap != u_map_alpha_Y.end()){
        if (lambda[mu] > 0){
            // how do we choose the population sizes, based on the lambdas?
            if (random_selection == 1){
                // let's sample the new sizes of the strains from a poissonian distribution with mean lambda[mu]
                boost::random::poisson_distribution<int> dist(lambda[mu]);
                imap->second = dist(RNG::rng);
            } else{
                // no fluctuations (n[mu] = lambda[mu])
                imap->second = int(lambda[mu]);
            }
            // go to the next entry
            imap++;
        } else{ // if lambda is zero, the strain is extinct
            // erase() will return the next iterator. We don't need to increment
            imap = u_map_alpha_Y.erase(imap);
        }
        mu++;
    }

    if (u_map_alpha_Y.size() == 0) {

        cout << "There isn't any individual left!\n";
        exit(EXIT_FAILURE);

    }

}

// -------------------- evolution: mutations -------------------- //
void mutation_step(u_map &u_map_alpha_Y, double s, double U_Y, double U_alpha, int R){

    // extract infos from u_map_alpha_Y
    vector<vector<double>> alpha_Y = get_alpha_Y(u_map_alpha_Y);
    vector<int> n = get_n(u_map_alpha_Y);

    int k_Y, k_lose, k_gain;
    int m;
    boost::random::binomial_distribution<int> dist;
    for (size_t mu = 0; mu < alpha_Y.size(); mu++){

        // fitness mutations
        dist.param(boost::random::binomial_distribution<int>::param_type(n[mu], U_Y));
        k_Y = dist(RNG::rng); // number of fitness mutations
        if(k_Y > 0)
            update_umap(fitness_mutation(alpha_Y[mu], k_Y, s, R), u_map_alpha_Y);

        // strategy mutations
        dist.param(boost::random::binomial_distribution<int>::param_type(n[mu], U_alpha));
        m = dist(RNG::rng); // number of strategy mutations (we have to divide them into gain and lose mutations)

        // gain
        dist.param(boost::random::binomial_distribution<int>::param_type(m, (double) R / (R + accumulate(alpha_Y[mu].begin(), alpha_Y[mu].end() - 1, 0.))));
        k_gain = dist(RNG::rng);
        if(k_gain > 0)
            update_umap(gain_mutation(alpha_Y[mu], k_gain, R), u_map_alpha_Y);

        // lose
        k_lose = m - k_gain;
        if(k_lose > 0)
            update_umap(lose_mutation(alpha_Y[mu], k_lose, R), u_map_alpha_Y);

        // we subtract the number of individuals with mutations from the original strain
        u_map_alpha_Y[alpha_Y[mu]] -= (k_Y + k_gain + k_lose); 

        if(u_map_alpha_Y[alpha_Y[mu]] == 0)   // if a strain is extinct
            u_map_alpha_Y.erase(alpha_Y[mu]); // we remove it

    }

}

vector<vector<double>>fitness_mutation(vector<double> alpha_Y_mu, int k_Y, double s, int R){

    vector<vector<double>> alpha_Y_new(k_Y, vector<double>(R + 1));

    alpha_Y_mu[R] += s; // we increase the (log) fitness by s

    for (int mu = 0; mu < k_Y; mu++){
        alpha_Y_new[mu] = alpha_Y_mu;
    }

    return alpha_Y_new;

}

vector<vector<double>> lose_mutation(vector<double> alpha_Y_mu, int k_lose, int R){

    vector<vector<double>> alpha_Y_new(k_lose, vector<double>(R + 1));
    vector<int> non_zeros(R + 1);
    int nz_size = 0;
    
    // we find the positions of the non-zero elements ** in alpha **
    for(int i = 0; i < R; i++){
        if(alpha_Y_mu[i]){ // alpha_Y_mu[i] != 0
            non_zeros[nz_size] = i;
            nz_size++;
        }
    }

    if(nz_size > 1){ // a strain must consume at least one resource

        boost::random::uniform_int_distribution<int> dist(0, nz_size - 1);
        vector<double> alpha_Y_copy(R + 1);

        for (int mu = 0; mu < k_lose; mu++){
            
            alpha_Y_copy = alpha_Y_mu;
            // select one resource with I_{\mu, i} = 1 and put I_{\mu, i} = 0
            alpha_Y_copy[non_zeros[dist(RNG::rng)]] = 0;
            alpha_Y_new[mu] = alpha_Y_copy;

        }

    } else {

        for (int mu = 0; mu < k_lose; mu++)
            alpha_Y_new[mu] = alpha_Y_mu;

    }
    
    return alpha_Y_new;

}

vector<vector<double>> gain_mutation(vector<double> alpha_Y_mu, int k_gain, int R){

    vector<vector<double>> alpha_Y_new(k_gain, vector<double>(R + 1));
    int num_res = accumulate(alpha_Y_mu.begin(), alpha_Y_mu.end() - 1, 0.); // number of resources consumed by the strain
    
    if (num_res < R){ // this slightly improves the performances

        boost::random::uniform_int_distribution<int> dist(0, R - 1);
        vector<double> alpha_Y_copy(R + 1);

        for (int mu = 0; mu < k_gain; mu++){

            alpha_Y_copy = alpha_Y_mu;
            // select one resource and put I_{\mu, i} = 1
            alpha_Y_copy[dist(RNG::rng)] = 1.;
            alpha_Y_new[mu] = alpha_Y_copy;

        }

    } else { // num_res = R

        for (int mu = 0; mu < k_gain; mu++)
            alpha_Y_new[mu] = alpha_Y_mu;

    }

    return alpha_Y_new;

}

void update_umap(const vector<vector<double>> &alpha_Y_new, u_map &u_map_alpha_Y){
    
    for(size_t mu = 0; mu < alpha_Y_new.size(); mu++)
        u_map_alpha_Y[alpha_Y_new[mu]]++;

}

// -------------------- extract infos from unordered_map -------------------- //
vector<vector<double>> get_alpha_Y(u_map &u_map_alpha_Y){

    vector<vector<double>> alpha_Y(u_map_alpha_Y.size(), vector<double>());

    size_t i = 0;
    for(auto const& imap: u_map_alpha_Y){
        alpha_Y[i] = imap.first;
        i++;    
    }
    return alpha_Y;

}

vector<vector<double>> get_alpha(u_map &u_map_alpha_Y, int R){

    vector<vector<double>> alpha_Y = get_alpha_Y(u_map_alpha_Y);
    vector<vector<double>> alpha(alpha_Y.size(), vector<double>(R));
    
    for (size_t mu = 0; mu < alpha_Y.size(); mu++)
        copy(alpha_Y[mu].begin(), alpha_Y[mu].end() - 1, alpha[mu].begin());

    return alpha;

}

vector<double> get_Y(u_map &u_map_alpha_Y, int R){

    vector<vector<double>> alpha_Y = get_alpha_Y(u_map_alpha_Y);
    vector<double> Y(alpha_Y.size());

    for (size_t mu = 0; mu < alpha_Y.size(); mu++)
        Y[mu] = alpha_Y[mu][R];

    return Y;

}

vector<int> get_n(u_map &u_map_alpha_Y){

    vector<int> n(u_map_alpha_Y.size());

    int i = 0;
    for(auto const& imap: u_map_alpha_Y){
        n[i] = imap.second;
        i++;
    }
    return n;

}

// -------------------- compute Lyapunov function -------------------- //
double lyapunov_function(u_map &u_map_alpha_Y, double chi, int R, vector<double> beta){

    double L = 0;
    double sum;

    for(auto const& imap: u_map_alpha_Y)
        L += imap.second * (1 + chi * accumulate(imap.first.begin(), imap.first.end() - 1, 0.));

    for(int i = 0; i < R; i++){
        sum = 0;
        for(auto const& imap: u_map_alpha_Y)
            sum += imap.first[i] * imap.second * exp(imap.first[R]);

        L -= beta[i] * log(sum);

    }

    return L;

}

// -------------------- summary: write on files -------------------- //
int summary(u_map u_map_alpha_Y, int D, int R, int n_start, int n_gen, int generalist, int random_selection, int write_on_file, double NU_Y, double NU_alpha, double eps, double s, double chi, vector<double> beta){

    u_map u_map_count_species;
    vector<vector<double>> ALPHA = get_alpha(u_map_alpha_Y, R);
    vector<int> N = get_n(u_map_alpha_Y);
    vector<double> Y = get_Y(u_map_alpha_Y, R);

    int sum_n = 0;
    vector<double> gen(R, 1); // generalist ecotype

    // let's create an unordered map with the abundances of the alphas (NOT alpha_Ys !)    
    // AKA an unordered map of the ecotypes (NOT of the strains)
    u_map u_map_n;
    for (size_t i = 0; i < ALPHA.size(); i++){

        u_map_n[ALPHA[i]] += N[i];
        sum_n += N[i]; // total population

    } 

    // ----- start of the sequencing ----- //

    // let's create the "inverse" map, with the cumulative sum of n's as keys
    // (this is useful to sample individuals knowing the abundances of the ecotypes)
    map<int, vector<double>> map_n;
    vector<vector<double>> alpha = get_alpha_Y(u_map_n); // in this case we obtain alpha using get_alpha_Y (IT'S NOT AN ERROR!)
    vector<int> n = get_n(u_map_n);
    partial_sum(n.begin(), n.end(), n.begin());

    for (size_t i = 0; i < n.size(); i++)
        map_n[n[i]] = alpha[i];

    // let's generate D random numbers between 1 and \sum n_\mu
    vector<int> rands(sum_n, 1);
    partial_sum(rands.begin(), rands.end(), rands.begin()); // rands = {1, 2, 3, ... , sum_n}
    for(auto& element: rands) element -= 1;                 // rands = {0, 1, 2, ... , sum_n - 1}
    shuffle(rands.begin(), rands.end(), RNG::rng);

    map<int, vector<double>>::iterator idx;

    // let's count the occurences of the sampled ecotypes
    for(int i = 0; i < min(D, sum_n); i++){ // beware! D could be greater than sum_n!
        idx = map_n.upper_bound(rands[i] - 1);
        u_map_count_species[idx->second]++;
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

    // ----- end of the sequencing ----- //

    // ----- we open a file to write on it information about the system ----- //

    char filename[100];
    FILE * fp;

    // if write_on_file == 1 we write the final infos (at the end of the run)
    // if write_on_file == 0 we write infos every dt_gen generations (checkpoints)
    if (write_on_file == 1)
        sprintf(filename, "./Data/run_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e.txt", (double) R, (double) n_start, (double) n_gen, (double) generalist, (double) random_selection, eps, s, NU_Y, NU_alpha, chi);
    else
        sprintf(filename, "./Data/n_gen_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e.txt", (double) R, (double) n_start, (double) n_gen, (double) generalist, (double) random_selection, eps, s, NU_Y, NU_alpha, chi);

    fp = fopen(filename, "a");

    // we write beta on the file
    for(int i = 0; i < R; i++){
        if (i < R - 1) fprintf(fp, "%f,", beta[i]);
        else           fprintf(fp, "%f\n", beta[i]);
    }

    // compute the mean fitness of the individuals that eat a resource: <e^Y>_i
    vector<double> mean_eY_i(R, 0.);
    double norm;

    for(int i = 0; i < R; i++){

        norm = 0;

        for(auto const& imap: u_map_alpha_Y){

            mean_eY_i[i] += imap.second * imap.first[i] * exp(imap.first[R]);
            norm += imap.second * imap.first[i];

        }

        mean_eY_i[i] /= norm;
        if (i < R - 1) fprintf(fp, "%f,", mean_eY_i[i]);
        else           fprintf(fp, "%f\n", mean_eY_i[i]);

    }

    // compute the mean LOG fitness of the individuals that eat a resource: <e^Y>_i
    vector<double> mean_Y_i(R, 0.);
    
    for(int i = 0; i < R; i++){

        norm = 0;

        for(auto const& imap: u_map_alpha_Y){

            mean_Y_i[i] += imap.second * imap.first[i] * imap.first[R];
            norm += imap.second * imap.first[i];

        }

        mean_Y_i[i] /= norm;
        if (i < R - 1) fprintf(fp, "%f,", mean_Y_i[i]);
        else           fprintf(fp, "%f\n", mean_Y_i[i]);

    }

    // compute F_i
    double F_i;
    for (int i = 0; i < R; i++){
        F_i = 0;
        for(auto &imap: u_map_n){
            F_i += imap.first[i] * imap.second / sum_n;
        }
    
        if (i < R - 1) fprintf(fp, "%f,", F_i);
        else		   fprintf(fp, "%f\n", F_i);
        
    }

    // compute the mean fitness of every ecotype: e^Y_b
    u_map fitness_eco;
    for (size_t i = 0; i < ALPHA.size(); i++)
        fitness_eco[ALPHA[i]] += exp(Y[i]) * (double) N[i] / u_map_n[ALPHA[i]];

    size_t mu = 0;
    for(auto const& imap: fitness_eco) {
        if (mu < fitness_eco.size() - 1) fprintf(fp, "%f,", imap.second);
        else                             fprintf(fp, "%f\n", imap.second);
        mu++;
    }

    // compute the mean LOG fitness of every ecotype: Y_b
    u_map LOG_fitness_eco;
    for (size_t i = 0; i < ALPHA.size(); i++)
        LOG_fitness_eco[ALPHA[i]] += Y[i] * (double) N[i] / u_map_n[ALPHA[i]];

    mu = 0;
    for(auto const& imap: LOG_fitness_eco) {
        if (mu < LOG_fitness_eco.size() - 1) fprintf(fp, "%f,", imap.second);
        else                             fprintf(fp, "%f\n", imap.second);
        mu++;
    }

    // compute the genome size of every ecotype (\sum_i a_{\sigma, i})
    double genome_size;
    mu = 0;
    for(auto const& imap: u_map_n) {
        genome_size = 0;
        for(int i = 0; i < R; i++)
            genome_size += imap.first[i];

        if (mu < u_map_n.size() - 1) fprintf(fp, "%f,", genome_size);
        else                         fprintf(fp, "%f\n", genome_size);
        mu++;
    }

    // we write the number of individuals of every ecotype
    n = get_n(u_map_n);
    for(size_t mu = 0; mu < n.size(); mu++){
        if(n[mu] > 0){
            if (mu < n.size() - 1) fprintf(fp, "%d,", n[mu]);
            else                   fprintf(fp, "%d\n", n[mu]);
        }
    }

    // we write the hashing value of every ecotype (to be able to recognize them later!)
    mu = 0;
    for(auto const& imap: u_map_n) {
        if (mu < u_map_n.size() - 1) fprintf(fp, "%li,", boost::hash_range(imap.first.begin(), imap.first.end()));
        else                         fprintf(fp, "%li\n", boost::hash_range(imap.first.begin(), imap.first.end()));
        mu++;
    }

    // compute the mean fitness of the entire population
    double mean_eY = 0; // <e^Y>
    for(auto const& imap: u_map_alpha_Y)
        mean_eY += exp(imap.first[R]) * imap.second / sum_n;

    fprintf(fp, "%f\n", mean_eY);

    // compute the mean LOG fitness of the entire population
    double mean_Y = 0; // <Y>
    for(auto const& imap: u_map_alpha_Y)
        mean_Y += imap.first[R] * imap.second / sum_n;

    fprintf(fp, "%f\n", mean_Y);

    // coefficient of variation CoV(e^Y) over the ecotypes
    double mean2 = 0;
    double mean = 0;
    
    for (size_t i = 0; i < alpha.size(); i++){ // for cycle over the ecotypes
        mean2 += fitness_eco[alpha[i]] * fitness_eco[alpha[i]] * u_map_n[alpha[i]] / sum_n;
        mean += fitness_eco[alpha[i]] * u_map_n[alpha[i]] / sum_n;
    }

    fprintf(fp, "%f\n", (mean2 - mean * mean) / (mean * mean));
    
    // coefficient of variation CoV(e^Y) over the strains
    mean2 = 0;
    mean = 0;
    for(auto const& imap: u_map_alpha_Y){

        mean2 += exp(imap.first[R]) * exp(imap.first[R]) * imap.second / sum_n;
        mean += exp(imap.first[R]) * imap.second / sum_n;

    }

    fprintf(fp, "%f\n", (mean2 - mean * mean) / (mean * mean));

    // coefficient of variation CoV(Y) over the ecotypes (LOG fitness)
    mean2 = 0;
    mean = 0;
    
    for (size_t i = 0; i < alpha.size(); i++){ // for cycle over the ecotypes
        mean2 += LOG_fitness_eco[alpha[i]] * LOG_fitness_eco[alpha[i]] * u_map_n[alpha[i]] / sum_n;
        mean += LOG_fitness_eco[alpha[i]] * u_map_n[alpha[i]] / sum_n;
    }

    fprintf(fp, "%f\n", (mean2 - mean * mean) / (mean * mean));

    // coefficient of variation CoV(Y) over the strains (LOG fitness)
    mean2 = 0;
    mean = 0;
    for(auto const& imap: u_map_alpha_Y) {

        mean2 += imap.first[R] * imap.first[R] * imap.second / sum_n;
        mean += imap.first[R] * imap.second / sum_n;

    }

    fprintf(fp, "%f\n", (mean2 - mean * mean) / (mean * mean));

    double k_mu = 0;
    for(auto const& imap: u_map_n)
        if(imap.first != gen) // = non generalist
            k_mu += 1. / inner_product(imap.first.begin(), imap.first.end(), imap.first.begin(), 0.) * imap.second;

    k_mu /= sum_n - u_map_n[gen];

    // if there's not a generalist ecotype, we have to remove it from the u_map
    if(u_map_n[gen] == 0)
        u_map_n.erase(gen);

    fprintf(fp, "%f,", u_map_n[gen] / sum_n); // generalist frequency
    fprintf(fp, "%f\n", R - k_mu); // R - <k_\mu>_{non generalist}

    fclose(fp);

    return u_map_count_species.size();

}

// -------------------- run the system -------------------- //
output_run run(int R, int n_start, int n_gen, int dt_gen, int D, double NU_Y, double NU_alpha, double eps, double s, int fast, int generalist, int random_selection, double chi){

    u_map u_map_alpha_Y;

    // initial population
    if (generalist == 1)
        initial_alpha_Y_gen(u_map_alpha_Y, R, n_start);

    else if (generalist == 0)
        initial_alpha_Y_spec(u_map_alpha_Y, R, n_start, s);

    else if (generalist == 2)
        complete_pool(u_map_alpha_Y, R, n_start, chi);

    else {cerr << "\nERROR: which initial population do you want?\n\n"; exit(EXIT_FAILURE);}
    
    vector<double> beta = compute_beta(R, eps); // resources' supply rates
    
    // probabilities of mutations
    double U_Y = NU_Y / n_start;
    double U_alpha = NU_alpha / n_start;

    int n_species;

    char filename[100];
    FILE * fp;

    // the core of the simulation is in this loop: mutations and selection!
    for (int i = 0; i < n_gen; i++){
        selection_step(u_map_alpha_Y, beta, R, n_start, random_selection, chi);
        mutation_step(u_map_alpha_Y, s, U_Y, U_alpha , R);

        if(i % dt_gen == 0 and fast == 0){ // checkpoints: saving information every dt_gen generations

            n_species = summary(u_map_alpha_Y, D, R, n_start, n_gen, generalist, random_selection, 0, NU_Y, NU_alpha, eps, s, chi, beta);

            // we open a file to write checkpoints on it
            sprintf(filename, "./Data/n_gen_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e.txt", (double) R, (double) n_start, (double) n_gen, (double) generalist, (double) random_selection, eps, s, NU_Y, NU_alpha, chi);
            fp = fopen(filename, "a");

            fprintf(fp, "%i\n", n_species);
            fprintf(fp, "%f\n", lyapunov_function(u_map_alpha_Y, chi, R, beta));

            fclose(fp);
         
        }

    }
    
    output_run result;
    result.scaled_rate = U_Y * s * s * (double) R / (U_alpha * eps * eps);
    result.n_species = summary(u_map_alpha_Y, D, R, n_start, n_gen, generalist, random_selection, 1, NU_Y, NU_alpha, eps, s, chi, beta);

    // let's compute and write on file the standard deviation of Y_i
    vector<double> Y_i = compute_eY_i(u_map_alpha_Y, beta, R);
    for (int i = 0; i < R; i++) Y_i[i] = log(Y_i[i]);
    double stdev = std_dev(Y_i);

    sprintf(filename, "./Data/run_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e.txt", (double) R, (double) n_start, (double) n_gen, (double) generalist, (double) random_selection, eps, s, NU_Y, NU_alpha, chi);
    fp = fopen(filename, "a");
    fprintf(fp, "%f\n", stdev);
    fclose(fp);
    
    return result;

}

vector<double> compute_eY_i(u_map &u_map_alpha_Y, vector<double> &beta, int R){

    vector<double> eY_i(R, 0);
    double small = 1E-11; // small number used to avoid divergence

    // let's compute \overline{Y}_i(t)
    for (int i = 0; i < R; i++){
        for(auto const& imap: u_map_alpha_Y)
            eY_i[i] += imap.first[i] * exp(imap.first[R]) * imap.second; // e^Y_i[i] += alpha[mu][i] * exp(Y[mu]) * n[mu]

        eY_i[i] /= beta[i];
        if(eY_i[i] == 0) eY_i[i] = small;
    }

    return eY_i;

}

double std_dev(vector<double> &v){

    double mean = accumulate(v.begin(), v.end(), 0.0) / v.size();

    double accum = 0.0;
    std::for_each (v.begin(), v.end(), [&](const double d) {
        accum += (d - mean) * (d - mean);
    });

    return sqrt(accum / (v.size() - 1)); // sample standard deviation

}

// -------------------- print data structures -------------------- //
void print_umap(u_map &umap){

    int last;

    for(auto &imap: umap){
        last = imap.first.size();
        for(int i = 0; i < last - 1; i++){
            cout << imap.first[i] << " ";
        }
        cout << exp(imap.first[last]);
        cout << " | " << imap.second << endl;

    }

}

void print_map(map<int, vector<double>> map){

    int last;

    for(auto &imap: map){

        last = imap.second.size();

        for(int i = 0; i < last - 1; i++){
            cout << imap.second[i] << " ";
        }

        cout << exp(imap.second[last]);
        cout << " | " << imap.first << endl;

    }

}

template <typename T>
void print_vec(vector<T> v){
    for(size_t i = 0; i < v.size(); i++) cout << v[i] << " ";
    cout << endl;
}