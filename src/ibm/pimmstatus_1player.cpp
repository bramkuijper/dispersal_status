// Gillespie simulation 
// for the evolution of age-dependent nongenetic effects

#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "bramauxiliary.h"

//#define NDEBUG

using namespace std;


// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

const size_t Npatch = 4000; 
const size_t Npp = 1;
const size_t numgen = 40000;
int sample = 20;
double d0 = 0.01; //some starting values
double h0 = 0.01;

int seed = -1; // the random seed
time_t total_time;  // track the total time the simulation is running
size_t generation = 0; // track the current generation (for debugging only)

size_t skip2 = numgen; // intervals of data output

double envt_switch[2] = {0,0}; // patch state switch rate
double k = 0.001; // cost of dispersal

// mortality of maladapted junior breeders in both envts
double c[2] = {0,0}; 

// mortality of adapted junior breeders in both envts
double C[2] = {0,0}; 

// mutation rates
double mu_h = 0.01;
double sdmu_h = 0.01;

// diploid individuals
struct Individual
{
    double d[2]; // one locus coding for migration
    double p1_imm[2]; // proportion z1 offspring by z1 immigrant
    double p2_imm[2]; // proportion z1 offspring by z2 immigrant
    double p1_nat[2]; // proportion z1 offspring by z1 native
    double p2_nat[2]; // proportion z1 offspring by z2 native
    bool z; // state of the individual wrt its environment, 0, 1

    bool immigrant; // individual has dispersed to the current patch
};

struct Patch
{
    Individual locals[Npp]; // the local senior
    bool state; //environmental state of the patch

};

Patch MetaPop[Npatch];


// statistics for sampling
// - number of patches in environment one or two
// - with individual in state (a,imm), (a,phil), (m, imm), (m, phil)
size_t Npatches[2][4];

// next to counts store lists with ids
size_t patch_ids[2][4][Npatch];

// give the outputfile a unique name
string filename("sim_disp");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// initialize the command line arguments to vary 
// the parameters
void init_arguments(int argc, char *argv[])
{
    envt_switch[0] = atof(argv[1]);
    envt_switch[1] = atof(argv[2]);
    c[0] = atof(argv[3]);
    c[1] = atof(argv[4]);
    C[0] = atof(argv[5]);
    C[1] = atof(argv[6]);
    k = atof(argv[7]);
    mu_h = atof(argv[8]);
    sdmu_h = atof(argv[9]);
    d0 = atof(argv[10]);
    h0 = atof(argv[11]);
}

// initialization function, runs at start of sim
// gives all the individuals some starting
// value
void init_pop()
{ 
    // start the time
    total_time = time(NULL);

    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // set the seed
    gsl_rng_set(r, seed)

    // initialize patch state counters and set them to 0
    for (int envt_i = 0; envt_i < 2; ++envt_i)
    {
        for (int n_states = 0; n_states <= 3; ++n_states)
        {
            Npatches[envt_i][n_states] = 0;
        }
    }

    // initialize all the patches
    for (size_t i = 0; i < Npatch; ++i)
    {
        //patches are randomly put in env + or -
        MetaPop[i].state = gsl_rng_uniform(r) < 0.5; 

        n_adapted = 0;
    
        //initialize breeders
        for (size_t j = 0; j < Npp; ++j)
        {
            // randomly assign initial phenotypes
            MetaPop[i].locals[j].z = gsl_rng_uniform(r) < 0.5;
            MetaPop[i].locals[j].immigrant = 0;

            if (MetaPop[i].locals[j].z == MetaPop[i].state)
            {
                ++n_adapted;
            }
        
            // initialize allelic values
            for (int allele = 0; allele < 2; ++allele)
            {
                MetaPop[i].locals[j].d[allele] = d0;
                MetaPop[i].locals[j].p1_imm[allele] = h0;
                MetaPop[i].locals[j].p2_imm[allele] = h0;
                MetaPop[i].locals[j].p1_nat[allele] = h0;
                MetaPop[i].locals[j].p2_nat[allele] = h0;
            }
        }

        // increase the counter of patches of
        // (i) environmental state 'state'
        // (iii) containing n_adapted individuals
        //
        size_t Npatches_of_this_type = 
            Npatches[MetaPop[i].state][0][Npp - n_adapted][0][n_adapted];

        // append this patch to the corresponding stack of patch id's
        patch_ids[MetaPop[i].state][0][Npp - n_adapted][0][n_adapted][Npatches_of_this_type] = i;

        // update counter
        ++Npatches[MetaPop[i].state][0][Npp - n_adapted][0][n_adapted];
    }
}


double mut_h(double p)
{
    assert(p >= 0 && p <= 1);
    p+= gsl_rng_uniform(r) < mu_h ? gsl_ran_gaussian(r,sdmu_h) : 0;
    p = p < 0 ? 0 : p > 1 ? 1 : p;

    return(p);
}

// make a kid out of parent asexually
// add it either to local offspring pool
// or to pool of dispersers
void create_kid(Individual &mother, Individual &Kid)
{
    // if it is a phenotype 1 mother, 
    // we'll need to express p1
    // otherwise p2
    double expr_h = 0.5;


    if (mother.immigrant)
    {
        // mother's phenotype is one
        expr_h = mother.z == 0 ? 0.5 * (mother.p1_imm[0] + mother.p1_imm[1]) 
                                    :
                                    0.5 * (mother.p2_imm[0] + mother.p2_imm[1]);
    }
    else
    {
        expr_h = mother.z == 0 ? 0.5 * (mother.p1_nat[0] + mother.p1_nat[1]) 
                                    :
                                    0.5 * (mother.p2_nat[0] + mother.p2_nat[1]);
    }

    // offspring gets z1 or z2 phenotype
    Kid.z = gsl_rng_uniform(r) < expr_h ? 0 : 1;

    // inherit alleles
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        Kid.d[allele_i] = mut_h(mother.d[allele_i]);
        Kid.p1_imm[allele_i] = mut_h(mother.p1_imm[allele_i]);
        Kid.p2_imm[allele_i] = mut_h(mother.p2_imm[allele_i]);
        Kid.p1_nat[allele_i] = mut_h(mother.p1_nat[allele_i]);
        Kid.p2_nat[allele_i] = mut_h(mother.p2_nat[allele_i]);
    }
}

// do the stats on the data
void write_data()
{
    // get variance and means
    double p1_imm,p2_imm,p1_nat,p2_nat,d,vard,varp1_imm,varp2_imm,varp1_nat,varp2_nat;
    double meanp1_imm = 0;
    double meanp2_imm = 0;
    double meanp1_nat = 0;
    double meanp2_nat = 0;
    double ssp1_imm = 0;
    double ssp2_imm = 0;
    double ssp1_nat = 0;
    double ssp2_nat = 0;

    double meand = 0;
    double ssd = 0;

    // loop through patches and individuals
    // and get stats on genotypes and patch freqs.
    for (size_t i = 0; i < Npatch; ++i)
    {
        for (size_t j = 0; j < Npp; ++j)
        {
            d = 0.5 * (MetaPop[i].locals[j].d[0] + MetaPop[i].locals[j].d[1]);
            p1_imm = 0.5 * (MetaPop[i].locals[j].p1_imm[0] + MetaPop[i].locals[j].p1_imm[1]);
            p2_imm = .5 * (MetaPop[i].locals[j].p2_imm[0] + MetaPop[i].locals[j].p2_imm[1]);
            p1_nat = 0.5 * (MetaPop[i].locals[j].p1_nat[0] + MetaPop[i].locals[j].p1_nat[1]);
            p2_nat = .5 * (MetaPop[i].locals[j].p2_nat[0] + MetaPop[i].locals[j].p2_nat[1]);

            meanp1_imm+=p1_imm;
            meanp2_imm+=p2_imm;
            meanp1_nat+=p1_nat;
            meanp2_nat+=p2_nat;
            meand+=d;
            ssp1_imm+=p1_imm*p1_imm;
            ssp2_imm+=p2_imm*p2_imm;
            ssp1_nat+=p1_nat*p1_nat;
            ssp2_nat+=p2_nat*p2_nat;
            ssd+=d*d;
        }    
    }

    DataFile << generation << ";"; 

    // write down all the patch frequency combinations
    for (size_t envt_i = 0; envt_i < 2; ++envt_i)
    {
        for (size_t n_adapted_imm = 0; n_adapted_imm <= Npp; ++n_adapted_imm)
        {
            for (size_t n_maladapted_imm = 0; n_maladapted_imm <= Npp - n_adapted_imm; ++n_maladapted_imm)
            {
                for (size_t n_maladapted_nat = 0; n_maladapted_nat <= Npp - n_adapted_imm - n_maladapted_imm; ++n_maladapted_nat)
                {
                    DataFile << setprecision(5) << (double)Npatches[envt_i][n_maladapted_imm][n_adapted_imm][n_maladapted_nat][Npp - n_maladapted_imm - n_adapted_imm - n_maladapted_nat] / Npatch << ";";
                }
            }
        }
    }

    meanp1_imm /= Npatch * 2;
    meanp2_imm /= Npatch * 2;
    varp1_imm = ssp1_imm / (Npatch * 2) - meanp1_imm*meanp1_imm;
    varp2_imm = ssp2_imm / (Npatch * 2) - meanp2_imm*meanp2_imm;
    meanp1_nat /= Npatch * 2;
    meanp2_nat /= Npatch * 2;
    varp1_nat = ssp1_nat / (Npatch * 2) - meanp1_nat*meanp1_nat;
    varp2_nat = ssp2_nat / (Npatch * 2) - meanp2_nat*meanp2_nat;
    meand /= Npatch * 2;
    vard = ssd / (Npatch * 2) - meand*meand;

    DataFile 
        << setprecision(5) << meanp1_imm << ";"
        << setprecision(5) << meanp2_imm << ";"
        << setprecision(5) << varp1_imm << ";"
        << setprecision(5) << varp2_imm << ";"
        << setprecision(5) << meanp1_nat << ";"
        << setprecision(5) << meanp2_nat << ";"
        << setprecision(5) << varp1_nat << ";"
        << setprecision(5) << varp2_nat << ";"
        << setprecision(5) << meand << ";"
        << setprecision(5) << vard << ";" << endl;
}

// headers of the datafile
void write_data_headers()
{
    DataFile << "generation;";

    for (size_t envt_i = 0; envt_i < 2; ++envt_i)
    {
        for (size_t n_adapted_imm = 0; n_adapted_imm <= Npp; ++n_adapted_imm)
        {
            for (size_t n_maladapted_imm = 0; n_maladapted_imm <= Npp - n_adapted_imm; ++n_maladapted_imm)
            {
                for (size_t n_maladapted_nat = 0; n_maladapted_nat <= Npp - n_adapted_imm - n_maladapted_imm; ++n_maladapted_nat)
                {
                    DataFile << setprecision(5) << "f_" << n_maladapted_imm << "_" << n_adapted_imm << "_" << n_maladapted_nat << "_" << n_adapted_nat << ";";
                }
            }
        }
    }
    
    DataFile << "meanp1_imm;meanp2_imm;varp1_imm;varp2_imm;meanp1_nat;meanp2_nat;varp1_nat;varp2_nat;meand;vard;" << endl;
}

// parameters at the end of the sim
void write_parameters()
{
    total_time = time(NULL) - total_time;
    DataFile << endl << endl << "d0;" << d0 << endl
                << "type;" << "phenotype-dependent" << endl
                << "p0;" << h0 << endl
                << "mu_h;" << mu_h << endl
                << "sdmu_h;" << sdmu_h << endl
                << "numgen;" << numgen << endl
                << "k;" << k << endl
                << "s1;" << envt_switch[0] << endl
                << "s2;" << envt_switch[1] << endl
                << "c1;" << c[0] << endl
                << "c2;" << c[1] << endl
                << "C1;" << C[0] << endl
                << "C2;" << C[1] << endl
                << "seed;" << seed << endl
                << "Npatches;" << Npatches << endl
                << "runtime;" << total_time << endl;
}


// switch the state of a patch in state
// patch_envt (e1 or e2)
// 0 <= n_adapted_imm <= Npp
// 0 <= n_maladapted_imm < Npp - n_adapted_imm
void switch_patch_state(size_t patch_envt, size_t n_adapted_imm, size_t n_maladapted_imm, size_t n_adapted_nat)
{
    // sample a random patch that fulfills the conditions
    size_t random_patch = 
        gsl_rng_uniform_int(r, 
            Npatches[patch_envt][n_maladapted_imm][n_adapted_imm][n_maladapted_nat][Npp - n_maladapted_imm - n_adapted_imm - n_maladapted_nat];
        );

    // get the patch_id 
    size_t random_patch_id = 
            patch_ids[patch_envt][n_maladapted_imm][n_adapted_imm][n_maladapted_nat][Npp - n_maladapted_imm - n_adapted_imm - n_maladapted_nat][random_patch];

    // error checking
    assert(MetaPop[random_patch_id].state == patch_envt);
   
#ifndef NDEBUG
    size_t n_adapted_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_adapted_check += MetaPop[random_patch_id].locals[breeder_i].z 
            == MetaPop[random_patch_id].state;
    }

    assert(n_adapted_check == n_adapted);
#endif

    // update patch characteristics
    bool patch_envt_new = !MetaPop[random_patch_id].state;
    MetaPop[random_patch_id].state = patch_envt_new;
        
        

    size_t n_adapted_new = Npp - n_adapted; 

    // update statistics
    //
    // delete patch id in the corresponding stack
    // by replacing it with the last patch id in the stack
    // and then deleting the last element (by reducing the counter by 1)
    patch_ids[patch_envt][n_adapted][random_patch] =
        patch_ids[patch_envt][n_adapted][
            --Npatches[patch_envt][n_adapted] 
        ];

    // add patch id to the correct stack
    patch_ids[patch_envt_new][n_adapted_new][
        Npatches[patch_envt_new][n_adapted_new]++
    ] = random_patch_id;
    
    // error checking
    assert(MetaPop[random_patch_id].state == !patch_envt);
   
#ifndef NDEBUG
    n_adapted_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_adapted_check += MetaPop[random_patch_id].locals[breeder_i].z 
            == MetaPop[random_patch_id].state;
    }

    assert(n_adapted_check == n_adapted_new);
#endif
}

// mortality of an individual in one of the patches in state 
// patch envt
//
// size_t patch_envt: the environment of the patch the mortality event took place in
// size_t n_adapted: a patch containing n_adapted adapted individuals
// bool mortality_maladapted: true, maladapted guy dies; false, adapted guy dies;

void mortality(size_t patch_envt, size_t n_adapted, bool mortality_maladapted)
{
    // sample a random patch that fulfills the conditions
    size_t random_patch = 
        gsl_rng_uniform_int(r, 
            Npatches[patch_envt][n_adapted]
        );

    // get the patch_id 
    size_t random_patch_id = 
        patch_ids[patch_envt][n_adapted][random_patch];

    // error checking
    assert(MetaPop[random_patch_id].state == patch_envt);
   
#ifndef NDEBUG
    size_t n_adapted_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_adapted_check += MetaPop[random_patch_id].locals[breeder_i].z 
            == MetaPop[random_patch_id].state;
    }

    assert(n_adapted_check == n_adapted);
#endif

    // fecundity selection absent
    // so whether we get immigrant 
    // only depends on dispersal probability
    
    size_t sampled = Npp + sample;

    // array to store in cumulative distribution over all sampled
    // individuals
    double choose_prob[sampled];

    // array to store patches from which the sampled individuals
    // come from
    size_t patch_origin[sampled];

    // array to store the ids of the individuals
    size_t individual_ids[sampled];

    // array with the individuals that may die
    size_t mortality_candidate_ids[Npp];
    size_t n_mortality_candidates = 0;


    double sumprob = 0;
    double d = 0;

    // make cumulative distribution of local competitive
    // ability over the locals (as dispersal is the only
    // trait affecting local competition just make a cumul
    // distribution of dispersal probabilities
    //
    // at the same time calculate which local breeders
    // might die
    for (size_t local_i = 0; local_i < Npp; ++local_i)
    {
        // calculate a local offspring's dispersal
        // capability
        d = 0.5 * (
                MetaPop[random_patch_id].locals[local_i].d[0] + 
                MetaPop[random_patch_id].locals[local_i].d[1]
                );

        // store the data
        patch_origin[local_i] = random_patch_id;
        individual_ids[local_i] = local_i;

        // store the cumulative distribution
        choose_prob[local_i] = sumprob + 1-d;
        sumprob = choose_prob[local_i];

        // store the mortality candidates
        if (
                (mortality_maladapted && 
                    MetaPop[random_patch_id].locals[local_i].z 
                        != MetaPop[random_patch_id].state
                )
                ||
                (!mortality_maladapted && 
                    MetaPop[random_patch_id].locals[local_i].z 
                        == MetaPop[random_patch_id].state
                )

            )
        {
            mortality_candidate_ids[n_mortality_candidates++] = local_i;
        }
    }

    // do some weighing of the dispersers
    double pi = 2.0 / sample;

    size_t id = 0;

    // make cumulative distribution over all sampled remote offspring
    // mimicking competition for breeding spots
    for (size_t remote_i = Npp; remote_i < sampled; ++remote_i)
    {
        // randomly pick patch
        patch_origin[remote_i] = gsl_rng_uniform_int(r,Npatch);

        // randomly pick individual within patch
        id = gsl_rng_uniform_int(r,Npp);

        individual_ids[remote_i] = id;

        // store individual id
        individual_ids[remote_i] = id;

        d =  0.5 * (
                        MetaPop[patch_origin[remote_i]].locals[id].d[0] + 
                        MetaPop[patch_origin[remote_i]].locals[id].d[1]
                        );

        choose_prob[remote_i] = sumprob + (1-k) * pi * d;
        sumprob = choose_prob[remote_i];
    }

    double deviate = gsl_rng_uniform(r) * sumprob;

    Individual Kid;

    bool kid_created = false;

    for (size_t i = 0; i < sampled; ++i)
    {
        if (deviate <= choose_prob[i])
        {
            assert(individual_ids[i] >= 0 && individual_ids[i] < Npp);
            assert(patch_origin[i] >= 0 && patch_origin[i] < Npatch);
            kid_created = true;

            create_kid(
                    MetaPop[patch_origin[i]].locals[individual_ids[i]],
                    Kid);
            break;
        }
    }

    assert(kid_created);

    // keep track of new state variables
    // to update counters (see below)
    size_t n_adapted_new = n_adapted;

    // then mortality 
    assert(n_mortality_candidates > 0 && n_mortality_candidates <= Npp);

    // sample a random candidate for mortality
    size_t candidate_id = 
        mortality_candidate_ids[gsl_rng_uniform_int(r, n_mortality_candidates)];

    assert(candidate_id >= 0 && candidate_id < Npp);

    assert(
            (!mortality_maladapted &&
            MetaPop[random_patch_id].locals[candidate_id].z ==
            MetaPop[random_patch_id].state
            )
            ||
            (mortality_maladapted &&
            MetaPop[random_patch_id].locals[candidate_id].z !=
            MetaPop[random_patch_id].state
            )
    );

    if (!mortality_maladapted)
    {
        --n_adapted_new;
    }

    // mortality and replacement
    MetaPop[random_patch_id].locals[candidate_id] = Kid;

    if (MetaPop[random_patch_id].locals[candidate_id].z ==
            MetaPop[random_patch_id].state)
    {
        ++n_adapted_new;
    }

    // update statistics
    //
    // delete patch id in the corresponding stack
    // by replacing it with the last patch id in the stack
    // and then deleting the last element (by reducing the counter by 1)
    patch_ids[patch_envt][n_adapted][random_patch] =
        patch_ids[patch_envt][n_adapted][
            --Npatches[patch_envt][n_adapted] 
        ];

    // add patch id to the correct stack
    patch_ids[patch_envt][n_adapted_new][
        Npatches[patch_envt][n_adapted_new]++
    ] = random_patch_id;

    // we're done.
    // error checking
    assert(MetaPop[random_patch_id].state == patch_envt);

#ifndef NDEBUG
    n_adapted_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_adapted_check += MetaPop[random_patch_id].locals[breeder_i].z 
            == MetaPop[random_patch_id].state;
    }

    assert(n_adapted_check == n_adapted_new);
#endif
}

// main function body
int main(int argc, char **argv)
{
    init_arguments(argc, argv);
    init_pop();
    write_parameters();
    write_data_headers();

    // there are this number of n probabilities:
    // 1. switching
    // 2. mortality 
    //
    // these can affect 2 * 2  different patch types
    // i. 2 x envt, ii. 3 x 0,1,2 adapted individuals per patch
    size_t nprobs = 2 * 2 * 3;

    double probs[nprobs]; // vector with probabilities

    bool done; // variable to get out of nested loop
    double prob; // actual event probability sampled from 
                // the cumulative prob distribution

    // set up as a markovian process of numgen^2 timesteps
    // if we would have used a single variable then we would 
    // run into problems as the minimum specification of the maximum size
    // of size_t = ~65500
    for (generation = 0; generation < numgen; ++generation)
    {
        for (size_t generation2 = 0; generation2 < 80000; ++generation2)
        {
            // generate cumulative prob. distribution of possible
            // events
            double sum_probs = 0;
            size_t prob_count = 0;

            // loop through the different patch types 
            for (size_t patch_envt = 0; patch_envt < 2; ++patch_envt)
            {
                // number of adapted individuals in a patch
                for (size_t n_adapted = 0; n_adapted <= 2; ++n_adapted)
                {
                    // probability that a patch switches state
                    sum_probs += envt_switch[patch_envt] 
                        * Npatches[patch_envt][n_adapted];

                    probs[prob_count++] = sum_probs; 

                    // mortality adapted
                    sum_probs += (double) n_adapted / Npp 
                        * 1.0 / (1.0 - C[patch_envt]) 
                        * Npatches[patch_envt][n_adapted];
                    
                    probs[prob_count++] = sum_probs; 

                    // mortality maladapted individual
                    sum_probs += (double) (Npp - n_adapted) / Npp 
                        * 1.0 / (1.0 - c[patch_envt]) 
                        * Npatches[patch_envt][n_adapted];
                    
                    probs[prob_count++] = sum_probs; 

                    assert(prob_count <= nprobs);
                }
            }

            // decide which event will happen
            //
            // draw value from cumul distribution
            prob = gsl_rng_uniform(r) *  sum_probs;

            prob_count = 0;

            // variable to break out of all the loops
            done = false;

            // locate the event
            // in the cumulative distribution
            for (size_t patch_envt = 0; patch_envt < 2; ++patch_envt)
            {
                if (done)
                {
                    break;
                }

                for (size_t n_adapted = 0; n_adapted <= 2; ++n_adapted)
                {
                    if (prob <= probs[prob_count])
                    {
                        switch_patch_state(patch_envt, n_adapted);
                        done = true;
                        break;
                    }

                    ++prob_count;

                    // mortality adapted
                    if (prob <= probs[prob_count])
                    {
                        mortality(patch_envt, n_adapted, false);
                        done = true;
                        break;
                    }

                    ++prob_count;

                    // mortality maladapted
                    if (prob <= probs[prob_count])
                    {
                        mortality(patch_envt, n_adapted, true);
                        done = true;
                        break;
                    }
                    
                    ++prob_count;
                }
            } // end for patch_envt
    
            if (generation2 % skip2 == 0)
            {
                write_data();
            }
        } //end for size_t generation2
    } // end for size_t _generation

    write_data();

    exit(1);
}
