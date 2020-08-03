#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <cmath>
#include <cassert>
#include "dispersal_status_2p.hpp"


// set a boundary for a value within (0,1)
double bound001(double val)
{
    return(val < 1e-07 ? 1e-07 : val > 1-(1e-07) ? 1-(1e-07) : val);
}

// set a boundary for a value as [0,1]
double bound01(double val)
{
    return(val < 0 ? 0 : val > 1 ? 1 : val);
}

// constructor
DispStatus2P::DispStatus2P(int argc, char **argv) :
    d{0.0}
    ,s{0.0,0.0}
    ,c{0.0,0.0}
    ,f{}
{
    // initialize the arguments
    init_arguments(argc, argv);

    // allocate a file to write data output to 
    std::ofstream output_file(filename.c_str());
    
    // write parameters to file first
    // after that the data
    // if we could write parameters last,
    // we would have a problem in case simulations
    // do typically not converge
    write_parameters(output_file);
    
    // write data headers
    write_data_headers(output_file);

    // current generation number
    long int timestep;

    // keep track whether everything has converged
    bool converged;

    // initialize the vector of 
    // traits in the next generation
    PatchFreq ftplus1{};

    double ptplus1[]{0.0,0.0,0.0,0.0};

    // now iterate the eco-evolutionary model
    for (timestep = 0; timestep < 1e07; ++timestep)
    {
        // first iterate the patch frequencies
        iterate_patch_frequencies(true);

        // then the relatedness coefficients (which rely on the 
        // right eigenvectors)
        iterate_relatedness();

        // write data every skip timesteps
        if (timestep % skip == 0)
        {
            write_data(output_file, timestep);
        }

        // we assume everything has converged
        // unless proven otherwise
        converged = true;
    }
} // DispStatus2P::DispStatus2P


// initialize the arguments from the command line
void DispStatus2P::init_arguments(int argc, char ** argv)
{
    s[0] = atof(argv[2]);
    s[1] = atof(argv[3]);
    c[0] = atof(argv[10]);
    c[1] = atof(argv[11]);
    d = atof(argv[11]);
}

// write values of the parameters
// to the output stream data_stream
void DispStatus2P::write_parameters(std::ofstream &data_stream)
{
    data_stream << std::endl << std::endl;
    data_stream << "s1;" << s[0] << std::endl;
    data_stream << "s2;" << s[1] << std::endl;
    data_stream << "d;" << d << std::endl;
    data_stream << "c1;" << c[0] << std::endl;
    data_stream << "c2;" << c[1] << std::endl;
}

// write headers to the datafile
void DispStatus2P::write_data_headers(std::ofstream &data_stream)
{
    data_stream << "time;";

    // translate the different states to a string
    std::string states[4] = {"d1","d2","n1","n2"};

    // headers for the patch frequencies 
    for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
    {
        for (int state1_idx = 0; state1_idx < 4; ++state1_idx)
        {
            for (int state2_idx = 0; state2_idx <= state1_idx; ++state2_idx)
            {
                data_stream << "f_e" << (envt_idx + 1) 
                    << "_" <<  states[state1_idx] 
                    << "_" <<  states[state2_idx] << ";";
            }
        }
    }

    data_stream << std::endl;
} // end DispStatus2P::write_data_headers()

// write the data to the output stream data_stream
void DispStatus2P::write_data(
        std::ofstream &data_stream, 
        size_t const time_val)
{
    data_stream  << time_val << ";";

    for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
    {
        for (int state1_idx = 0; state1_idx < 4; ++state1_idx)
        {
            for (int state2_idx = 0; state2_idx <= state1_idx; ++state2_idx)
            {
                data_stream << f(envt_idx, state1_idx, state2_idx) << ";";
            }
        }
    }


    data_stream << std::endl;
} // end DispStatus2P::write_data


// iterate relatedness until equilibrium 
// is achieved
void DispStatus2P::iterate_relatedness()
{
} // end iterate relatedness

// calculates the probability that an already dead individual in state 
// state_dead_individual will be replaced by an individual in exactly
// the same state (in which case no change in patch frequency occurs)
double DispStatus2P::replace_mort_with_same(
        State const state_dead_individual // state of the dead individual
        ,bool const envt_local_patch // environment in local patch
        ,State const state_other_individual)
{
    // individual who replaces dead individual
    // is philopatric, meaning either the dead individual
    // or the other individual gave birth to it
    if (state_dead_individual == n1)
    {
        return((1.0 - d) * (p[state_dead_individual] + p[state_other_individual]));
    }
    else if (state_dead_individual == n2)
    {
        return((1.0 - d) * (1.0 - p[state_dead_individual] + 1.0 - p[state_other_individual]));
    }

    // otherwise indivual who replaces dead individual
    // is immigrant and we have to calculate overall probability
    // by summing over all patches
    double prob_replace = 0.0;

    for (int envt_remote_idx = 0; envt_remote_idx < 2; ++envt_remote_idx)
    {
        for (int state1_remote_idx = 0; state1_remote_idx < 4; ++state1_remote_idx)
        {
            for (int state2_remote_idx = 0; state2_remote_idx < state1_remote_idx; 
                    ++state2_remote_idx)
            {
                prob_replace += d * (
                        state_dead_individual == d1 ?
                            p[state1_remote_idx] + p[state2_remote_idx] // z1
                            :
                            1.0 - p[state1_remote_idx] + 1.0 - p[state2_remote_idx]); // z2
            }
        }
    }

    assert(prob_replace >= 0.0);
    assert(prob_replace <= 1.0);

    return(prob_replace);
} // end  DispStatus2P::replace_mort_with_same()


double DispStatus2P::mort(bool const envt, State const state)
{
    if (envt)
    {
        return(state == n2 || state == d2 ? 1.0 : c[envt]);
    }

    return(state == n1 || state == d1 ? 1.0 : c[envt]);
} // DispStatus2P::mort()


void DispStatus2P::mort_influx(
        bool const envt_state_target
        ,State const state1_target
        ,State const state2_target)
{
    State state1_idx_s,state2_idx_s;

    // first thing that can happen, state1_target individual dies
    // and state2_individual == state2_target
    // then just calculate recruitment of state1 individual == state1_target
    //
    // if however, state2_individual == state1_target
    // then we need to recruit individual that matches state2_target
    
    double total_rate = 0.0;
        
    for (int state1_idx = 0; state1_idx < 4; ++state1_idx)
    {
        state1_idx_s = static_cast<State>(state1_idx);

        for (int state2_idx = 0; state2_idx <= state1_idx; ++state2_idx)
        {
            state2_idx_s = static_cast<State>(state2_idx);

            if (state2_idx_s == state2_target && state1_idx_s == state1_target)
            {
                continue;
            }

            if (state2_idx_s == state2_target)
            {
                total_rate += mort(state1_idx_s) * 
                    f(envt_state_target, state1_idx_s, state2_idx_s)

                    // TODO
        }
    }

    double total_rate = mort(

} // end DispStatus2P::mort_influx()

// iterate patch frequencies until convergence
void DispStatus2P::iterate_patch_frequencies(bool output)
{
    // auxiliary variable to see whether stuff has been converged
    bool converged;

    // initialize the patch frequencies for the next timestep
    PatchFreq ftplus1{};

    State state1_idx_s, state2_idx_s;

    // now iterate the patch frequencies, reproductive values
    // and relatedness coefficients until convergence
    for (size_t iter_time_step = 0; 
            iter_time_step < tmax_ecological;
            ++iter_time_step)
    {
        for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
        {
            for (int state1_idx = 0; state1_idx < 4; ++state1_idx)
            {
                state1_idx_s = static_cast<State>(state1_idx);


                for (int state2_idx = 0; state2_idx <= state1_idx; ++state2_idx)
                {
                    state2_idx_s = static_cast<State>(state2_idx);

                    ftplus1(envt_idx, state1_idx_s, state2_idx_s) = 

                        // outflux because environment changes away from envt_idx
                        -s[envt_idx] * f(envt_idx, state1_idx_s, state2_idx_s)

                        // influx because patches in opposite environmental state
                        // change to envt_idx
                        + s[1 - envt_idx] * f(1 - envt_idx, state1_idx_s, state2_idx_s)
                        
                        // outflux because of mortality events of player 1
                        -mort(envt_idx, state1_idx_s) * f(envt_idx, state1_idx_s, state2_idx_s) 
                        * (1.0 - replace_mort_with_same(state1_idx_s, envt_idx, state2_idx_s))

                        // outflux because of mortality events of player 2
                        -mort(envt_idx, state2_idx_s) * f(envt_idx, state1_idx_s, state2_idx_s)
                        * (1.0 - replace_mort_with_same(state2_idx_s, envt_idx, state1_idx_s));
                    
                        // influx because of mortality events elsewhere
                        + mort_influx(envt_idx, state1_idx_s, state2_idx_s);
                }
            }
        }
    }// end for loop

}

