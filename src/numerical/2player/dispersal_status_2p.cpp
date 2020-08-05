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
    ,p{0.0,0.0,0.0,0.0}
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
    // patch frequencies in the next generation
    PatchFreq ftplus1{};

    // initialize vector of epigenetic switch rates in the next
    // generation
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
    s[0] = atof(argv[1]);
    s[1] = atof(argv[2]);
    c[0] = atof(argv[3]);
    c[1] = atof(argv[4]);
    d = atof(argv[5]);
    p[0] = atof(argv[6]);
    p[1] = atof(argv[7]);
    p[2] = atof(argv[8]);
    p[3] = atof(argv[9]);
    filename = argv[6];
}

// write values of the parameters
// to the output stream data_stream
void DispStatus2P::write_parameters(std::ofstream &data_stream)
{
    data_stream << "s1;" << s[0] << std::endl;
    data_stream << "s2;" << s[1] << std::endl;
    data_stream << "d;" << d << std::endl;
    data_stream << "c1;" << c[0] << std::endl;
    data_stream << "c2;" << c[1] << std::endl;
    data_stream << std::endl << std::endl;
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
        return((1.0 - d) * 0.5 * (p[state_dead_individual] + p[state_other_individual]));
    }
    else if (state_dead_individual == n2)
    {
        return((1.0 - d) * (1.0 - 0.5 *(p[state_dead_individual] + p[state_other_individual])));
    }

    // otherwise indivual who replaces dead individual
    // is immigrant and we have to calculate overall probability
    // by summing over all patches
    double prob_replace = 0.0;

    // all options for recruiting dead individual from remote site
    for (int envt_remote_idx = 0; envt_remote_idx < 2; ++envt_remote_idx)
    {
        for (int state1_remote_idx = 0; state1_remote_idx < 4; ++state1_remote_idx)
        {
            for (int state2_remote_idx = 0; state2_remote_idx <= state1_remote_idx; 
                    ++state2_remote_idx)
            {
                prob_replace += d * f(envt_remote_idx, state1_remote_idx, state2_remote_idx) * (
                        state_dead_individual == d1 ?
                            0.5 * (p[state1_remote_idx] + p[state2_remote_idx]) // z1
                            :
                            1.0 - 0.5 *(p[state1_remote_idx] + p[state2_remote_idx])); // z2
            }
        }
    }

    std::cout << "dead: " << state_dead_individual << "," << envt_local_patch << "," << state_other_individual << ")" << prob_replace << std::endl;

    assert(prob_replace >= 0.0);
    assert(prob_replace <= 1.0);

    return(prob_replace);
} // end  DispStatus2P::replace_mort_with_same()

// recruit an immigrant from a remote patch who is in state state_target
double DispStatus2P::recruit_immigrant(State const state_target)
{
    // auxiliary variables containing the current states
    // of breeders in a remote patch (envt_idx, state1_idx_s, state2_idx_s)
    State state1_idx_s;
    State state2_idx_s;

    // state_target should always be d1 or d2 as individual
    // should be an immigrant (we are recruiting from a remote patch
    // after all)
    assert(state_target == d1 || state_target == d2);

    // now calculate total probability of recruiting such 
    // an individual from a remote patch
    double total_prob = 0.0;
    for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
    {
        for (int state1_idx = 0; state1_idx < 4; ++state1_idx)
        {
            state1_idx_s = static_cast<State>(state1_idx);

            for (int state2_idx = 0; state2_idx <= state1_idx; ++state2_idx)
            {
                state2_idx_s = static_cast<State>(state2_idx);

                total_prob += d * f(envt_idx, state1_idx_s, state2_idx_s) *
                    (state_target == d1 ? 0.5 * (p[state1_idx_s] + p[state2_idx_s]) : 
                     1.0 - 0.5 *(p[state1_idx_s] + p[state2_idx_s]));
            }
        }
    }

    return(total_prob);
} // end DispStatus2P::recruit_immigrant()

double DispStatus2P::mort(bool const envt, State const state)
{
    if (envt)
    {
        return(state == n2 || state == d2 ? 1.0 : c[envt]);
    }

    return(state == n1 || state == d1 ? 1.0 : c[envt]);
} // DispStatus2P::mort()


// calculates rate of influx due to mortalities
// to patches in state envt_state_target
// with breeders state1_target, state2_target
double DispStatus2P::mort_influx(
        bool const envt_state_target
        ,State const state1_target
        ,State const state2_target)
{
    State state1_idx_s,state2_idx_s;
    double p_recruit_state1_target, p_recruit_state2_target;

    // first thing that can happen, state1 individual dies
    // and state2_individual == state2_target
    // then calculate recruitment of newly recruited individual == state1_target
    //
    // if however, state2_individual == state1_target
    // then we need to recruit individual that matches state2_target

    // Repeat the same for the case that state2 individual dies
    double total_rate = 0.0;
       
    for (int state1_idx = 0; state1_idx < 4; ++state1_idx)
    {
        state1_idx_s = static_cast<State>(state1_idx);

        for (int state2_idx = 0; state2_idx <= state1_idx; ++state2_idx)
        {
            state2_idx_s = static_cast<State>(state2_idx);

            // if state of origin is the same as target state
            // there is no effective influx/outflux so we can ignore this option
            if (state2_idx_s == state2_target && state1_idx_s == state1_target)
            {
                continue;
            }

            // preliminary work 1:
            // calculate probability  that a state1_target
            // newborn will be recruited
            p_recruit_state1_target = 0.0;
            
            // replacement is a dispersing individual
            if (state1_target == d1 || state1_target == d2)
            {
                p_recruit_state1_target = recruit_immigrant(state1_target);
            }
            else // replacement is a philopatric individual
            {
                p_recruit_state1_target = (1.0 - d) * (state1_target == n1 ? 
                        0.5 * (p[state1_idx_s] + p[state2_idx_s])
                        :
                        1.0 - 0.5 * (p[state1_idx_s] + p[state2_idx_s]));
            }
            
            // preliminary work 2:
            // calculate probability that a state2_target
            // newborn will be recruited
            p_recruit_state2_target = 0.0;

            // replacement is a dispersing individual
            if (state2_target == d1 || state2_target == d2)
            {
                p_recruit_state2_target = recruit_immigrant(state2_target);
            }
            else // replacement is a philopatric individual
            {
                p_recruit_state2_target = (1.0 - d) * (state2_target == n1 ? 
                        0.5 * (p[state1_idx_s] + p[state2_idx_s])
                        :
                        1.0 - 0.5 *(p[state1_idx_s] + p[state2_idx_s]));
            }

            // individual in state1 dies
            // calculate probabilities of replacement

            // if state2_idx == state2_target then existing
            // state2 breeder will remain state2 breeder after
            // mortality. Hence, we 
            // need to recruit a new individual
            // in state state1_target
            if (state2_idx_s == state2_target)
            {
                // prob individual with state1 dies, 
                // x probability of recruitment
                total_rate += mort(envt_state_target, state1_idx_s) * 
                    f(envt_state_target, state1_idx_s, state2_idx_s)
                    * p_recruit_state1_target;
            }
            
            // if state2_idx == state1_target then 
            // state2 breeder will become the new state1
            // breeder after mortality. 
            //
            // We then
            // need to recruit a state2 individual
            // in state state2_target
            // note that state2_idx == state2_target AND state2_idx == state1_target
            // can both be true at the same time
            // hence no else if but just an if in the statement below
            if (state2_idx_s == state1_target)
            {
                total_rate += mort(envt_state_target, state1_idx_s) * 
                    f(envt_state_target, state1_idx_s, state2_idx_s)
                    * p_recruit_state2_target;
            }

            // individual in state2 dies
            //
            //
            // if state1_idx == state1_target then existing
            // state1 breeder will remain state1 breeder after
            // mortality. Hence, we 
            // need to recruit a new individual
            // in state state2_target
            if (state1_idx_s == state1_target)
            {
                total_rate += mort(envt_state_target, state2_idx_s) *
                    f(envt_state_target, state1_idx_s, state2_idx_s)
                    * p_recruit_state2_target;
            }

            // if state1_idx == state2_target then existing
            // state1 breeder will become state2 breeder after
            // mortality. Hence, we 
            // need to recruit a new individual
            // in state state1_target
            if (state1_idx_s == state2_target)
            {
                total_rate += mort(envt_state_target, state2_idx_s) *
                    f(envt_state_target, state1_idx_s, state2_idx_s)
                    * p_recruit_state1_target;
            }
        } // end for state2_idx
    } // end for state1_idx;

    return(total_rate);
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
        std::cout << "time: " << iter_time_step << std::endl;
        for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
        {
            for (int state1_idx = 0; state1_idx < 4; ++state1_idx)
            {
                state1_idx_s = static_cast<State>(state1_idx);

                for (int state2_idx = 0; state2_idx <= state1_idx; ++state2_idx)
                {
                    state2_idx_s = static_cast<State>(state2_idx);

                    ftplus1(envt_idx, state1_idx_s, state2_idx_s) = 
                        f(envt_idx, state1_idx_s, state2_idx_s) + eul * (

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
                        * (1.0 - replace_mort_with_same(state2_idx_s, envt_idx, state1_idx_s))
                    
                        // influx because of mortality events elsewhere
                        + mort_influx(envt_idx, state1_idx_s, state2_idx_s)

                        // and we are done.
                            );

                    if (ftplus1(envt_idx, state1_idx_s, state2_idx) != 0.0)
                    {
                        assert(std::isnormal(ftplus1(envt_idx, state1_idx_s, state2_idx)));
                    }
                }
            }
        }

        // assume convergence by default
        converged = true;

        // now update patch frequency values and check for convergence
        for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
        {
            for (int state1_idx = 0; state1_idx < 4; ++state1_idx)
            {
                state1_idx_s = static_cast<State>(state1_idx);

                for (int state2_idx = 0; state2_idx <= state1_idx; ++state2_idx)
                {
                    state2_idx_s = static_cast<State>(state2_idx);

                    if (
                            fabs(ftplus1(envt_idx, state1_idx_s, state2_idx_s) - 
                                f(envt_idx, state1_idx_s, state2_idx_s)) 
                            > ecology_vanish_bound
                                )
                    {
                        converged = false;
                    }

                    std::cout << "f(" << envt_idx << "," << state1_idx_s << "," << state2_idx_s << ") " 
                        << ftplus1(envt_idx, state1_idx_s, state2_idx_s) << " "  
                        << f(envt_idx, state1_idx_s, state2_idx_s) << std::endl;

                    f(envt_idx, state1_idx_s, state2_idx_s) = ftplus1(envt_idx, state1_idx_s, state2_idx_s);

                }
            }
        }

        if (converged)
        {
            break;
        }

    }// end for loop iter_time_step

} // end DispStatus2P::iterate_patch_frequencies()

