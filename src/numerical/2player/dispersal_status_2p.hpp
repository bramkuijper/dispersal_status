#ifndef DISP_STATUS_2P_HPP_ 
#define DISP_STATUS_2P_HPP_

#include <vector>
#include "patchfreq.hpp"


class DispStatus2P 
{
    private:
        
        // parameters
        double d; // dispersal rate
        double s[2]; // environmental switching

        // the evolving switch rates
        double p[4]; 

        static constexpr int tmax_ecological = 1e07;
        static constexpr double ecology_vanish_bound = 1e-08;
        static constexpr double eul = 1e-02;
        // cost of local adaptation in envt 0, 1
        double c[2];

        // vector to store the patch frequencies
        // there are 10 possible patch frequencies
        PatchFreq f;
        
        // store filename for output
        std::string filename;


        // skip data output until every skipth generation
        size_t skip = 10;
        
        void init_arguments(int argc, char **argv);

        // iterate the eigenvectors until convergence
        void iterate_patch_frequencies(bool output);

        // iterate relatdness coefficients until convergence
        void iterate_relatedness();

        // write headers to the data file
        void write_data_headers(std::ofstream &data_stream);

        // write data to the data file
        void write_data(std::ofstream &data_stream, size_t const generation_i);
        // write parameters to the data file
        void write_parameters(std::ofstream &data_stream);

        // recruit an immigrant in state state
        double recruit_immigrant(State const state);

        // mortality of an individual in state state
        // breeding in environment envt
        double mort(bool const envt, State const state);

        // calculates rate of influx due to mortalities
        // to patches in state envt_state_target
        // with breeders state1_target, state2_target
        double mort_influx(
            bool const envt_state_target
            ,State const state1_target
            ,State const state2_target);
        
        // calculates the probability that an already dead individual in state 
        // state_dead_individual will be replaced by an individual in exactly
        // the same state (in which case no change in patch frequency occurs)
        double replace_mort_with_same(
                State const state_dead_individual // state of the dead individual
                ,bool const envt_local_patch // environment in local patch
                ,State const state_other_individual);
    public:

        // class constructor
        DispStatus2P(int argc, char **argv);
};

#endif
