#ifndef DISP_STATUS_2P_HPP_ 
#define DISP_STATUS_2P_HPP_

#include <vector>


enum State
{
    Female = 0,
    Male = 1
};




class DispStatus2P 
{
    private:
        // parameters
        double d; // dispersal rate
        double sigma[2]; // environmental switching

        // cost of local adaptation in envt 0, 1
        double c[2];

        // vector to store the patch frequencies
        // there are 10 possible patch frequencies
        std::vector <double> f;
        
        // reproductive values 
        std::vector <double> v;
        
        std::vector <double> v;

        // store filename for output
        std::string filename;

        // for any ecological variable x in (u,v,Q)
        // we assume convergence when |x_t+1 - x_t| < 
        // ecology_vanish_bound
        double ecology_vanish_bound = 1e-08;

        // skip data output until every skipth generation
        size_t skip = 10;
        
        // resident fitness
        double w_resident(
                bool const envt_j
                ,Sex const sex_j
                ,bool const envt_i
                ,Sex const sex_i
                );

        double F(
                bool const envt_i);
        
        double dFdsflocff(
                bool const envt_deriv
                ,bool const envt_i);

        double dFdsmlocmm(
                bool const envt_deriv
                ,bool const envt_i);
        
        double dFdqremote(
                int const combn_i_deriv
                ,bool const envt_i);

        double dFdqloc(
                int const combn_i_deriv
                ,bool const envt_i);

        // probability of producing offspring with
        // phenotype offspring_phenotype
        double px(
                bool const offspring_phenotype
                ,bool const envt_female
                ,bool const envt_male);

        double dpxdsf(
                bool const envt_deriv // sf[envt_deriv] is the signal over which we take derivs
                ,bool const offspring_phenotype
                ,bool const envt_female
                ,bool const envt_male);

        double dpxdsm(
                bool const envt_deriv // sf[envt_deriv] is the signal over which we take derivs
                ,bool const offspring_phenotype
                ,bool const envt_female
                ,bool const envt_male);

        double dpxdq(
                int const q_deriv // q[q_deriv] is the offspring responsiveness over which we take derivs
                ,bool const offspring_phenotype
                ,bool const envt_female
                ,bool const envt_male);


        // total number of competiting juveniles
        // of sex_j
        // in a patch that transitions
        // from envt_i to envt_j,
        double C(
                bool const envt_i // origin envt_i
                ,bool const envt_j // future envt_j
                ,Sex const sex_j);
        
        // derivative of the 
        // total number of competiting juveniles
        // of sex_j
        // in a patch that transitions
        // from envt_i to envt_j,
        //
        // with respect to sf[envt_deriv]
        double dCdsflocff(
                bool const envt_deriv
                ,bool const envt_i // origin envt_i
                ,bool const envt_j // future envt_j
                ,Sex const sex_j); // futur sex_j

        double dCdsmlocmm(
                bool const envt_deriv
                ,bool const envt_i // origin envt_i
                ,bool const envt_j // future envt_j
                ,Sex const sex_j); // futur sex_j
        
        double dCdqremote(
                int const combn_i_deriv // which q to take derivatives of
                ,bool const envt_i // origin envt_i
                ,bool const envt_j // future envt_j
                ,Sex const sex_j); // futur sex_j

        double dCdqlocal(
                int const combn_i_deriv // which q to take derivatives of
                ,bool const envt_i // origin envt_i
                ,bool const envt_j // future envt_j
                ,Sex const sex_j); // futur sex_j

        double envt_switch(
                bool const envt_i
                ,bool const envt_j);

        // juvenile survival 
        double omega(
                bool const offspring_phenotype
                ,bool const envt_j);

        // global frequency of patches in environmental state envt_j
        double f(bool const envt_j);

        // probability nonlocally sired individual stays at the local patch
        double hn(
                    Sex const sex
                    , bool const envt_i
                    , bool const envt_j);
        
        // probability locally sired individual stays at the local patch
        double hl(
                    Sex const sex
                    , bool const envt_i
                    , bool const envt_j);

        double QJlocal_local(bool const envt_i);
        double QJnonlocal_local(bool const envt_i);
        double QJnonlocal_nonlocal(bool const envt_i);

        // coefficient of consanguinity in adults
        // in the next generation
        double Qtplus1(
                Sex const sex1
                ,Sex const sex2
                ,bool const envt_i);


        // selection gradient on female signal
        double dWdsf(bool const envt);
        
        // derivative of the transition matrix B
        // element b_(dim_i, dim_j) wrt sflocff[envt_deriv]
        double dWdsflocff(
                int dim_i
                ,int dim_j
                ,bool envt_deriv);
        
        
        // derivative of the transition matrix B
        // element b_(dim_i, dim_j) wrt sflocmf[envt_deriv]
        double dWdsflocmf(
                int dim_i
                ,int dim_j
                ,bool envt_deriv);

        // derivative of the transition matrix B
        // element b_(dim_i, dim_j) wrt sffocenvt_deriv]
        double dWdsffoc(
                int dim_i
                ,int dim_j
                ,bool envt_deriv);

        
        // selection gradient on male signal
        double dWdsm(bool const envt);
        
        // derivative of the transition matrix B
        // element b_(dim_i, dim_j) wrt smlocmm[envt_deriv]
        double dWdsmlocmm(
                int dim_i
                ,int dim_j
                ,bool envt_deriv);
        
        
        // derivative of the transition matrix B
        // element b_(dim_i, dim_j) wrt dWdsmlocfm[envt_deriv]
        double dWdsmlocfm(
                int const dim_i
                ,int const dim_j
                ,bool const envt_deriv);

        // derivative of the transition matrix B
        // element b_(dim_i, dim_j) wrt smfoc]
        double dWdsmfoc(
                int const dim_i
                ,int const dim_j
                ,bool const envt_deriv);

        double dWdq(int const combn_i);
        
        double dWdqfocf(
                int const dim_i // future
                ,int const dim_j // current
                ,int const combn_i_deriv // which of the 4 responsiveness traits to take derivs over
                );
        
        double dWdqfocm(
                int const dim_i // future
                ,int const dim_j // current
                ,int const combn_i_deriv // which of the 4 responsiveness traits to take derivs over
                );

        double dWdqlocf(
                int const dim_i // future
                ,int const dim_j // current
                ,int const combn_i_deriv // which of the 4 responsiveness traits to take derivs over
                );
        
        double dWdqlocm(
                int const dim_i // future
                ,int const dim_j // current
                ,int const combn_i_deriv // which of the 4 responsiveness traits to take derivs over
                );
        
        double dWdqfocfremote(
                int const dim_i // future
                ,int const dim_j // current
                ,int const combn_i_deriv // which of the 4 responsiveness traits to take derivs over
                );
        
        double dWdqfocmremote(
                int const dim_i // future
                ,int const dim_j // current
                ,int const combn_i_deriv // which of the 4 responsiveness traits to take derivs over
                );
        
        double dWdqlocfremote(
                int const dim_i // future
                ,int const dim_j // current
                ,int const combn_i_deriv // which of the 4 responsiveness traits to take derivs over
                );
        
        double dWdqlocmremote(
                int const dim_i // future
                ,int const dim_j // current
                ,int const combn_i_deriv // which of the 4 responsiveness traits to take derivs over
                );

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

    public:

        // class constructor
        MatPat(int argc, char **argv);
};

#endif
