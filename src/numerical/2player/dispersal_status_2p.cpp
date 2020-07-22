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
DispStatus2P::DispStatus2P(int argc, char **argv):
    d{0.0}
    ,sigma{0.0,0.0}
    ,c{0.0,0.0}
    ,f(10, 1.0/10)
    ,v(20, 1.0/10)

{
    // initialize the arguments
    init_arguments(argc, argv);

    std::ofstream output_file(filename.c_str());
    
    // write parameters to file first
    // if convergence takes long, we can already summarize
    // over all output files and see how well we are doing
    write_parameters(output_file);
    
    // write data headers
    write_data_headers(output_file);

    // current generation number
    long int timestep;

    // keep track whether everything has converged
    bool converged;

    // initialize the vector of 
    // traits in the next generation
    double ptplus1[]{0.0,0.0,0.0,0.0};
    double sftplus1[]{0.0,0.0};
    double qtplus1[]{0.0,0.0,0.0,0.0};

    // now iterate the eco-evolutionary model
    for (timestep = 0; timestep < 1e07; ++timestep)
    {
        // first iterate the eigenvectors 
        iterate_patch_frequencies(true);

        // then the relatedness coefficients (which rely on the 
        // right eigenvectors)
        iterate_relatedness();

        // write data every skip timesteps
        if (generation % skip == 0)
        {
            write_data(output_file, generation);
        }

        // we assume everything has converged
        // unless proven otherwise
        converged = true;

        // now assess the paternal signals in both envts 
        for (size_t envt_i = 0; envt_i < 2; ++envt_i)
        {
            // signal value in next timestep is bounded,
            // and the sum of the previous value, plus the 
            // selection gradient times genetic variance
            smtplus1[envt_i] = bound001(sm[envt_i] + Vsm * dWdsm(envt_i));
            sftplus1[envt_i] = bound001(sf[envt_i] + Vsf * dWdsf(envt_i));

            // check for numerical errors (NaN's)
            if (std::isnan(smtplus1[envt_i]) != 0)
            {
                std::cout << "range error in sm" << envt_i << std::endl;
                assert(std::isnan(smtplus1[envt_i]) == 0);
            }
            
            // check for numerical errors (NaN's)
            if (std::isnan(sftplus1[envt_i]) != 0)
            {
                std::cout << "range error in sf" << envt_i << std::endl;
                assert(std::isnan(sftplus1[envt_i]) == 0);
            }

            // check whether values have converged
            // this can be because there is little difference
            // between hi_t+1 and hi_t
            // or because hi_t never evolves beyond 0.02 after
            // a substantial amount of time
            if (!(fabs(smtplus1[envt_i] - sm[envt_i]) < 1e-08 
                    ||
                (generation > 1e04 && (sm[envt_i] < 0.001 || sm[envt_i] > 0.999))))
            {
                converged = false;
            }
            
            if (!(fabs(sftplus1[envt_i] - sf[envt_i]) < 1e-08 
                    ||
                (generation > 1e04 && (sf[envt_i] < 0.001 || sf[envt_i] > 0.999))))
            {
                converged = false;
            }
        } // end for size envt_i

        // now values of offspring responsiveness
        for (size_t combn_i = 0; combn_i < 3; ++combn_i)
        {
            // helping value in next timestep is bounded,
            // and the sum of the previous value, plus the 
            // selection gradient times genetic variance
            qtplus1[combn_i] = bound001(q[combn_i] + Vq * dWdq(combn_i));

            // check for numerical errors (NaN's)
            if (std::isnan(qtplus1[combn_i]) != 0)
            {
                std::cout << "range error in q" << combn_i << std::endl;
                assert(std::isnan(qtplus1[combn_i]) == 0);
            }
            
            // check whether values have converged
            // this can be because there is little difference
            // between hi_t+1 and hi_t
            // or because hi_t never evolves beyond 0.02 after
            // a substantial amount of time
            if (!(fabs(qtplus1[combn_i] - q[combn_i]) < 1e-08 
                    ||
                (generation > 1e04 && (q[combn_i] < 0.001 || q[combn_i] > 0.999))))
            {
                converged = false;
            }
            
        }

        // all values of h_i have converged, end the iteration
        if (converged)
        {
            break;
        }
        
        // no convergence, update values
        // and continue to the next step of the iteration
        for (int envt_i = 0; envt_i < 2; ++envt_i)
        {
            sf[envt_i] = sftplus1[envt_i];
            sm[envt_i] = smtplus1[envt_i];
        }
//        
        for (int combn_i = 0; combn_i < 3; ++combn_i)
        {
            q[combn_i] = qtplus1[combn_i];
        }
    
    }// end for generaiton i 

    // write the last bit of data
    if (generation % skip != 0)
    {
        write_data(output_file, generation);
    }
    
    // we are done
    output_file.close();

    // in case you want to explore code further, start at selection gradients
    // below, i.e., dWdh
}


// initialize the arguments from the command line
void MatPat::init_arguments(int argc, char ** argv)
{
    l = atof(argv[1]);
    sigma[0] = atof(argv[2]);
    sigma[1] = atof(argv[3]);
    n[Female] = atof(argv[4]);
    n[Male] = atof(argv[5]);
    d[Female] = atof(argv[6]);
    d[Male] = atof(argv[7]);
    error[Female] = atof(argv[8]);
    error[Male] = atof(argv[9]);
    c[0] = atof(argv[10]);
    c[1] = atof(argv[11]);
    y[0] = atof(argv[12]);
    y[1] = atof(argv[13]);
    sf[0] = atof(argv[14]);
    sf[1] = atof(argv[15]);
    sm[0] = atof(argv[16]);
    sm[1] = atof(argv[17]);
    q[0] = atof(argv[18]);
    q[1] = atof(argv[19]);
    q[2] = atof(argv[20]);
    Vsm = atof(argv[21]);
    Vsf = atof(argv[22]);
    Vq = atof(argv[23]);

    for (int sex_i = 0; sex_i < 2; ++sex_i)
    {
        for (int sex_j = 0; sex_j < 2; ++sex_j)
        {
            for (int envt_i = 0; envt_i < 2; ++envt_i)
            {
                Q[sex_i][sex_j][envt_i] = 0.1;
            }
        }
    }
}

void MatPat::write_parameters(std::ofstream &data_stream)
{
    data_stream << std::endl << std::endl;
    data_stream << "l;" << l << std::endl;
    data_stream << "sigma1;" << sigma[0] << std::endl;
    data_stream << "sigma2;" << sigma[1] << std::endl;
    data_stream << "nm;" << n[Male] << std::endl;
    data_stream << "nf;" << n[Female] << std::endl;
    data_stream << "dm;" << d[Male] << std::endl;
    data_stream << "df;" << d[Female] << std::endl;
    data_stream << "Em;" << error[Male] << std::endl;
    data_stream << "Ef;" << error[Female] << std::endl;
    data_stream << "Vsm;" << Vsm << std::endl;
    data_stream << "Vsf;" << Vsf << std::endl;
    data_stream << "Vq;" << Vq << std::endl;
    data_stream << "c1;" << c[0] << std::endl;
    data_stream << "c2;" << c[1] << std::endl;
    data_stream << "y1;" << y[0] << std::endl;
    data_stream << "y2;" << y[1] << std::endl;

    for (int i = 0; i < 3; ++i)
    {
        data_stream << "q" << i << "_init"  << q[i] << std::endl;
    }

    for (int i = 0; i < 2; ++i)
    {
        data_stream << "sf" << i << "_init" << sf[i] << std::endl;
    }
    
    for (int i = 0; i < 2; ++i)
    {
        data_stream << "sm" << i << "_init" << sm[i] << std::endl;
    }
}

// write headers to the datafile
void MatPat::write_data_headers(std::ofstream &data_stream)
{
    data_stream << "generation;";
   
    // headers for the stable class frequencies 
    // headers for the reproductive values
    for (size_t ev_i = 0; ev_i < u.size(); ++ev_i)
    {
        data_stream << "u" << ev_i << ";" << "v" << ev_i << ";";
    }

    // headers for the eigenvalues
    data_stream << "ev" << ";";

    for (size_t sex_i = 0; sex_i < 2; ++sex_i)
    {
        for (size_t sex_j = 0; sex_j < 2; ++sex_j)
        {
            for (size_t envt_i = 0; envt_i < 2; ++envt_i)
            {
                data_stream << "Q" 
                    << (sex_i == Female ? "f" : "m") 
                    << (sex_j == Female ? "f" : "m")
                    << envt_i << ";";
            }
        }
    }

    for (size_t envt_i = 0; envt_i < 2; ++envt_i)
    {
        data_stream << "sf" << envt_i << ";" << "sm" << envt_i << ";";
    }
    
    for (size_t combn_i = 0; combn_i < 3; ++combn_i)
    {
        data_stream << "q" << combn_i << ";";
    }

    data_stream << std::endl;
}


void MatPat::write_data(std::ofstream &data_stream, size_t const generation_i)
{
    data_stream  << generation_i << ";";

    for (size_t ev_i = 0; ev_i < u.size(); ++ev_i)
    {
        data_stream << u[ev_i] << ";" << v[ev_i] << ";";
    }

    data_stream << ev << ";";

    for (int sex_i = 0; sex_i < 2; ++sex_i)
    {
        for (int sex_j = 0; sex_j < 2; ++sex_j)
        {
            for (int envt_i = 0; envt_i < 2; ++envt_i)
            {
                data_stream << 
                    std::setprecision(8) << Q[static_cast<Sex>(sex_i)][static_cast<Sex>(sex_j)][envt_i] << ";";
            }
        }
    }

    for (size_t envt_i = 0; envt_i < 2; ++envt_i)
    {
        data_stream << sf[envt_i] << ";" << sm[envt_i] << ";";
    }
    
    for (size_t combn_i = 0; combn_i < 3; ++combn_i)
    {
        data_stream << q[combn_i] << ";";
    }

    data_stream << std::endl;
}


double MatPat::QJlocal_local(bool const envt_i)
{
    double Qval = 
        0.25 * (1.0 / n[Female] + (n[Female] - 1.0)/n[Female] 
                * Q[Female][Female][envt_i])
        + 0.5 * Q[Female][Male][envt_i]
        + 0.25 * (1.0 / n[Male] + (n[Male] - 1.0)/n[Male] 
                * Q[Male][Male][envt_i]);

    return(Qval);
}

double MatPat::QJnonlocal_local(bool const envt_i)
{
    double Qval = 
        0.25 * (1.0 / n[Female] + (n[Female] - 1.0)/n[Female] 
                * Q[Female][Female][envt_i])
        + 0.25 * Q[Female][Male][envt_i];

    return(Qval);
}



double MatPat::QJnonlocal_nonlocal(bool const envt_i)
{
    double Qval = 
        0.25 * (1.0 / n[Female] + 
                (n[Female] - 1.0)/n[Female] 
                * Q[Female][Female][envt_i]);
    
    return(Qval);
}


// probability of samping a philopatric offspring 
// sired by a nonlocal male
double MatPat::hn(
        Sex const sex
        , bool const envt_i
        , bool const envt_j)
{
    double hn_value = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        double sum_remote_male_probabilities = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_remote_male_probabilities += 
                f(envt_remote) * px(z_off_i, envt_i, envt_remote);
        }

        hn_value += n[Female]* F(envt_i) * (1.0 - l)
            * sum_remote_male_probabilities
            * (1.0 - d[sex]) * omega(z_off_i, envt_j) / 
               C(envt_i, envt_j, sex);
    }

    return(hn_value);
}

// probability of samping a philopatric offspring 
// sired by a local male
double MatPat::hl(
        Sex const sex
        , bool const envt_i
        , bool const envt_j)
{
    double hl_value = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        hl_value += n[Female]* F(envt_i) * l
            * px(z_off_i, envt_i, envt_i)
            * (1.0 - d[sex]) * omega(z_off_i, envt_j) / 
               C(envt_i, envt_j, sex);
    }

    return(hl_value);
}


// coefficient of consanguity between two adult breeders
// one of sex1, the other of sex2, both living in 
// environment envt_i
double MatPat::Qtplus1
(
        Sex const sex1
        ,Sex const sex2
        ,bool const envt_i
)
{
    double Qt1 = 0.0;

    // calculate all ways in which current environment is 1
    double total_prob_envt_is_i = 0.0;

    for (int envt_y = 0; envt_y < 2; ++envt_y)
    {
        total_prob_envt_is_i += f(envt_y) * envt_switch(envt_y, envt_i);
    }

    for (int envt_x = 0; envt_x < 2; ++envt_x)
    {
        Qt1 += f(envt_x) * envt_switch(envt_x, envt_i) / total_prob_envt_is_i

            // hl x hl
            * (
                    hl(sex1, envt_x, envt_i) * hl(sex2, envt_x, envt_i) * 
                            QJlocal_local(envt_x)

                    // hl x hn
                    + hl(sex1, envt_x, envt_i) * hn(sex2, envt_x, envt_i) * 
                            QJnonlocal_local(envt_x)
                    
                    // hl x hn
                    + hn(sex1, envt_x, envt_i) * hl(sex2, envt_x, envt_i) * 
                            QJnonlocal_local(envt_x)

                    // hn x hn
                    + hn(sex1, envt_x, envt_i) * hn(sex2, envt_x, envt_i) *
                            QJnonlocal_nonlocal(envt_x)
                );
    }

    assert(std::isinf(Qt1) == 0);
    assert(std::isnan(Qt1) == 0);

    return(Qt1);                    
}

// iterate relatedness until equilibrium 
// is achieved
void MatPat::iterate_relatedness()
{
    bool converged;

    // vector to store t+1 values
    double Qtplus1vec[2][2][2];

    // now iterate Q until convergence
    for (size_t iter_time_step = 0; 
            iter_time_step < 1e07;
            ++iter_time_step)
    {
        converged = true;

        for (int sex_i = 0; sex_i < 2; ++sex_i)
        {
            for (int sex_j = 0; sex_j < 2; ++sex_j)
            {
                for (int envt_i = 0; envt_i < 2; ++envt_i)
                {
                    Sex _sex_i = static_cast<Sex>(sex_i);
                    Sex _sex_j = static_cast<Sex>(sex_j);
//
                    // update Qt+1
                    Qtplus1vec[sex_i][sex_j][envt_i] = 
                        bound01(Qtplus1(_sex_i, _sex_j, envt_i));

                    // check for convergence
                    if (fabs(Qtplus1vec[sex_i][sex_j][envt_i] - Q[sex_i][sex_j][envt_i]) 
                            > ecology_vanish_bound)
                    {
                        converged = false;
                    }
                }
            }
        }

        if (converged)
        {
            break;
        }
        
        // update relatedness values
        for (int sex_i = 0; sex_i < 2; ++sex_i)
        {
            for (int sex_j = 0; sex_j < 2; ++sex_j)
            {
                for (int envt_i = 0; envt_i < 2; ++envt_i)
                {
                    Q[sex_i][sex_j][envt_i] = Qtplus1vec[sex_i][sex_j][envt_i];
                }
            }
        }
    } // end iterate
}

// iterate the system of eigenvectors until convergence
// we need to iterate the system as the eigenvector itself
// is a function of the patch frequencies, which are the right eigenvec
void DispStatus2P::iterate_patch_frequencies(bool output)
{
    // dimensions are envt_i * male = 2^3 = 8
    size_t ndim = 4;

    // allocate GSL things to work on the eigenvectors
    //
    // allocate a vector to calculate eigenvalues
    gsl_vector_complex *eval = gsl_vector_complex_alloc(ndim);

    // allocate a matrix to calculate eigenvectors
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(ndim, ndim);
     
    // allocate a vector to calculate eigenvalues for the tranpose
    gsl_vector_complex *evalT = gsl_vector_complex_alloc(ndim);
    
    // allocate a matrix to calculate eigenvectors for the tranpose
    gsl_matrix_complex *evecT = gsl_matrix_complex_alloc(ndim, ndim);

    gsl_matrix *m = gsl_matrix_alloc(ndim,ndim);
    gsl_matrix *mT = gsl_matrix_alloc(ndim,ndim);
    gsl_matrix *m_stats = gsl_matrix_alloc(ndim,ndim);
     
    gsl_eigen_nonsymmv_workspace * workspace_w = gsl_eigen_nonsymmv_alloc(ndim);
    gsl_eigen_nonsymmv_workspace * workspace_wT = gsl_eigen_nonsymmv_alloc (ndim);

    // auxiliary variable to see whether stuff has been converged
    bool converged;

    // auxiliary variables to store the lhs
    double v_i_tplus1, u_i_tplus1, ev_tplus1;

    double Atmp;

    size_t row_i, col_j;
        
    // now iterate the eigenvectors
    for (size_t iter_time_step = 0; 
            iter_time_step < 1e07;
            ++iter_time_step)
    {
        row_i = 0;
        col_j = 0;

        //std::cout << iter_time_step << std::endl;
        // generate the entries of the resident matrix
        for (size_t sex_i = 0; sex_i < 2; ++sex_i)
        {
            for (size_t envt_i = 0; envt_i < 2; ++envt_i)
            {
                for (size_t sex_j = 0; sex_j < 2; ++sex_j)
                {
                    for (size_t envt_j = 0; envt_j < 2; ++envt_j)
                    {
                        Atmp = w_resident(
                                envt_i
                                ,static_cast<Sex>(sex_i)
                                ,envt_j
                                ,static_cast<Sex>(sex_j));

                        gsl_matrix_set(m, row_i, col_j, Atmp);
                        gsl_matrix_set(mT, col_j, row_i, Atmp);
                        gsl_matrix_set(m_stats, row_i, col_j, Atmp);

                        ++col_j;

                        if (col_j > 3)
                        {
                            col_j = 0;
                        }
                        
                    }// end for envt_j
                } // end for sex_j

                ++row_i;

                if (row_i > 3)
                {
                    row_i = 0;
                }
            } // end for envt_i
        } //end for sex_i
        
        // now calculate evs
        // set up the workspaces
        gsl_eigen_nonsymmv(m, eval, evec, workspace_w);
        gsl_eigen_nonsymmv(mT, evalT, evecT, workspace_wT);
                
        

        gsl_eigen_nonsymmv_sort (eval, evec, 
                                GSL_EIGEN_SORT_ABS_DESC);
            
        // get the first column from the eigenvector matrix
        // giving us the right eigenvector
        gsl_vector_complex_view evec_i 
                   = gsl_matrix_complex_column(evec, 0);

        // get the first column from the eigenvector matrix
        // of the transpose of the left eigenvector
        // giving us the left eigenvector
        gsl_eigen_nonsymmv_sort (evalT, evecT, 
                                GSL_EIGEN_SORT_ABS_DESC);

        // get the last column from the eigenvector
        gsl_vector_complex_view evecT_i 
                   = gsl_matrix_complex_column(evecT, 0);

        double product_vu = 0;
        double total_u = 0;

        // make variables that normalize both eigenvectors
        for (size_t iter1 = 0; iter1 < ndim; ++iter1)
        {
            // norm for right ev
            total_u += fabs(GSL_REAL(
                   gsl_vector_complex_get(&evec_i.vector,iter1)));

            // norm for left ev
            product_vu += fabs(GSL_REAL(
                   gsl_vector_complex_get(&evec_i.vector,iter1)))
                   * fabs(GSL_REAL(
                   gsl_vector_complex_get(&evecT_i.vector,iter1)));
        }


        product_vu /= total_u;
       
        // by default we assume convergence until proved otherwise
        converged = true;
        
        // get the results from the eigenvector
        for (size_t iter1 = 0; iter1 < ndim; ++iter1)
        {
            // get the right eigenvalue out of the eigenvector of A
            u_i_tplus1 = fabs(GSL_REAL(
                        gsl_vector_complex_get(&evec_i.vector,iter1)))/
                total_u;
            
            // get the left eigenvalue out of the eigenvector of A^T
            v_i_tplus1 = fabs(GSL_REAL(
                        gsl_vector_complex_get(&evecT_i.vector,iter1)))/
                product_vu;

            if (fabs(u[iter1] - u_i_tplus1) > ecology_vanish_bound)
            {
                converged = false;
            }

            if (fabs(v[iter1] - v_i_tplus1) > ecology_vanish_bound)
            {
                converged = false;
            }

            v[iter1] = v_i_tplus1;
            u[iter1] = u_i_tplus1;

            assert(std::isnan(v[iter1]) == 0);
            assert(std::isnan(u[iter1]) == 0);

        }
        
        // now get the eigenvalue
        ev_tplus1 = GSL_REAL(gsl_vector_complex_get(eval, 0));

        if (fabs(ev - ev_tplus1) > ecology_vanish_bound)
        {
            converged = false;
        }

        ev = ev_tplus1;

        if (converged)
        {
            if (output)
            {
//                std::string filename_ev = filename + "_ev";
//                std::ofstream output_file_ev(filename_ev.c_str());
//
//                for (size_t ev_i = 0; ev_i < ndim; ++ev_i)
//                {
//                    for (size_t ev_j = 0; ev_j < ndim; ++ev_j)
//                    {
//                        output_file_ev << "w" << ev_i << "_" << ev_j << ";";
//                    }
//                }
//                
//                output_file_ev << std::endl;
//
//                for (size_t ev_i = 0; ev_i < ndim; ++ev_i)
//                {
//                    for (size_t ev_j = 0; ev_j < ndim; ++ev_j)
//                    {
//                        output_file_ev << gsl_matrix_get(m_stats, ev_i, ev_j) << ";";
//                    }
//                }
//
//                output_file_ev << std::endl;
//
//                output_file_ev.close();

            }

            break;
        }
        
//        for (size_t iter = 0; iter < u.size(); ++iter)
//        {
//            std::cout << "u" << iter << ": " << u[iter] << std::endl;
//        }
//        
//        for (size_t iter = 0; iter < u.size(); ++iter)
//        {
//            std::cout << "v" << iter << ": " << v[iter] << std::endl;
//        }
//        std::cout << "ev:" << ev << std::endl;
    } // end iteration


    // free all the stuff we've used
    gsl_eigen_nonsymmv_free (workspace_w);
    gsl_eigen_nonsymmv_free (workspace_wT);
    gsl_vector_complex_free(evalT);
    gsl_matrix_complex_free(evecT);
    gsl_matrix_free(mT);
    gsl_matrix_free(m);
    gsl_matrix_free(m_stats);

}

//// the selection gradient on helping behaviour
//double MatPat::dWdh(size_t const k)
//{
//    double dwdhk = 0;
//    double rloc, dwfoc, dwloc;
//
//    // sum over all rows of the mutant
//    // transition matrix
//    for (size_t i = 0; i <= nhmax; ++i)
//    {
//        // sum over all columns of the 
//        // mutant transition matrix
//        for (size_t j = 0; j <= nhmax; ++j)
//        {
//            // calculate relatedness of mother to
//            // locally born offspring
//            rloc = 1.0 / n + (n - 1.0) / n * Q[j];
//
//            // calculate selection on the focal mother
//            dwfoc = dwdahfoc(i, j, k);
//            // calculate selection on any local mother
//            dwloc = dwdahloc(i, j, k);
//
//            // check for numerical errors
//            assert(std::isnan(rloc) == 0);
//            assert(std::isnan(dwfoc) == 0);
//            assert(std::isnan(dwloc) == 0);
//
//            // We now have v_i * u_j * D[b_ij,h_k]
//            // where bij is the entry of the mutant transition
//            // matrix
//            dwdhk += v[i] * u[j] * (dwfoc + rloc * dwloc);
//        }
//    }
//
//    return(dwdhk);
//}

// get entries of the resident transition matrix
double MatPat::w_resident(
        bool const envt_j // future environment
        ,Sex const sex_j // future sex
        ,bool const envt_i // current environment
        ,Sex const sex_i // current sex
        )
{
    double w = 0.0;

    // first calculate f
    if (sex_i == Female)
    {
        // loop through all possible offspring phenotypes that could
        // be produced
        for (int z_off = 0; z_off < 2; ++z_off)
        {
            // w through nonlocal mating
            double w_nonlocal_mating = 0.0;

            // w through dispersal
            double w_dispersal = 0.0;

            // calculate payoffs through nonlocal mating
            // and dispersal
            for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
            {
                w_nonlocal_mating += f(envt_remote) * px(z_off, envt_i, envt_remote);

                w_dispersal += f(envt_remote) * envt_switch(envt_remote, envt_j) * 
                    n[sex_j] * d[sex_j] / 
                    C(envt_remote, envt_j, sex_j);
            }

            w += F(envt_i) * 
                (
                    l * px(z_off, envt_i, envt_i)
                    + (1.0 - l) * w_nonlocal_mating
                ) * omega(z_off, envt_j)
                * (envt_switch(envt_i, envt_j) * n[sex_j] * 
                 (1.0 - d[sex_j]) / C(envt_i, envt_j, sex_j)
                + w_dispersal);
        }

        return(0.5 * w);
    }

    for (int z_off = 0; z_off < 2; ++z_off)
    {
        // calculate remotely sired remote offspring
        double remote_remote_offspring = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            remote_remote_offspring += 
                f(envt_remote) * envt_switch(envt_remote, envt_j) * d[sex_j] * n[sex_j]
                / C(envt_remote, envt_j, sex_j);
        }

        double remote_offspring_local_mating = 0.0;
        double remote_offspring_nonlocal_mating = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            remote_offspring_local_mating += 
                f(envt_remote) * n[sex_j] * d[sex_j] * envt_switch(envt_remote, envt_j) /
                C(envt_remote, envt_j, sex_j);

            remote_offspring_nonlocal_mating += 
                f(envt_remote) * px(z_off, envt_remote, envt_i)
                            * F(envt_remote) * omega(z_off, envt_j)
                            * (envt_switch(envt_remote, envt_j) * n[sex_j] * (1.0 - d[sex_j]) / 
                                    C(envt_remote, envt_j, sex_j)
                                + remote_remote_offspring);
        }

        w += F(envt_i) * l * n[Female] / n[Male] *
            (
                 px(z_off, envt_i, envt_i) * omega(z_off, envt_j)
                 * (
                     n[sex_j] * (1.0 - d[sex_j]) / 
                        C(envt_i, envt_j, sex_j)
                         * envt_switch(envt_i, envt_j)
                     + remote_offspring_local_mating))
            + n[Female] / n[Male] * (1.0 - l) * remote_offspring_nonlocal_mating;
    }

    return(0.5 * w);
}

// probability of producing offspring with
// phenotype offspring_phenotype
double MatPat::px(
    bool const offspring_phenotype
    ,bool const envt_female
    ,bool const envt_male)
{
    double psignal_mom = (1.0 - error[Female]) * sf[envt_female]
        + error[Female] * sf[!envt_female];

    double psignal_dad = (1.0 - error[Male]) * sm[envt_male]
         + error[Male] * sm[!envt_male];

    assert(psignal_mom >= 0.0);
    assert(psignal_mom <= 1.0);
    assert(psignal_dad >= 0.0);
    assert(psignal_dad <= 1.0);

    double prob_offspring_z1 = psignal_mom * psignal_dad * q[2]
        + (
                psignal_mom * (1.0 - psignal_dad) 
                + 
                (1.0 - psignal_mom) * psignal_dad
        ) * q[1]
        + (1.0 - psignal_mom) * (1.0 - psignal_dad) * q[0];

    assert(prob_offspring_z1 >= 0.0);
    assert(prob_offspring_z1 <= 1.0);

    // prob of producing offspring phenotype z2 is requested
    if (offspring_phenotype)
    {
        return(1.0 - prob_offspring_z1);
    }

    return(prob_offspring_z1);
} // end double MatPat::px(

double MatPat::dCdqremote(
        int combn_i_deriv 
        ,bool const envt_i
        ,bool const envt_j
        ,Sex const sex_j)
{
    double cval = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // information from a male who mates at the focal patch
        // but who is from a nonlocal patch himself
        double sum_male_nonlocal_mating = 0.0;
        double d_sum_male_nonlocal_mating_dq = 0.0;

        for (int remote_envt_i = 0; remote_envt_i < 2; ++remote_envt_i)
        {
            sum_male_nonlocal_mating += 
                f(envt_i) * px(z_off_i, envt_i, remote_envt_i);

            d_sum_male_nonlocal_mating_dq +=
                f(envt_i) * dpxdq(combn_i_deriv, z_off_i, envt_i, remote_envt_i);
        }

        cval += // philopatric contribution to competing juveniles
            n[Female] * (1.0 - d[sex_j]) * (
                dFdqremote(combn_i_deriv, envt_i) * 
                ( l * px(z_off_i, envt_i, envt_i)
                  + (1.0 - l) * sum_male_nonlocal_mating
                ) 
                +
                F(envt_i) *
                    (1.0 - l) * d_sum_male_nonlocal_mating_dq
            ) * omega(z_off_i, envt_j);
    }

    return(cval);
} // end double MatPat::dCdqremote

double MatPat::dCdqlocal(
        int combn_i_deriv 
        ,bool const envt_i
        ,bool const envt_j
        ,Sex const sex_j)
{
    double cval = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // information from a male who mates at the focal patch
        // but who is from a nonlocal patch himself
        double sum_male_nonlocal_mating = 0.0;

        for (int remote_envt_i = 0; remote_envt_i < 2; ++remote_envt_i)
        {
            sum_male_nonlocal_mating += 
                f(envt_i) * px(z_off_i, envt_i, remote_envt_i);
        }

        cval += // philopatric contribution to competing juveniles
            n[Female] * (1.0 - d[sex_j]) * (
                dFdqloc(combn_i_deriv, envt_i) * 
                ( l * px(z_off_i, envt_i, envt_i)
                  + (1.0 - l) * sum_male_nonlocal_mating
                ) 
                +
                F(envt_i) *
                    l * dpxdq(combn_i_deriv, z_off_i, envt_i, envt_i)
            ) * omega(z_off_i, envt_j);
    }

    return(cval);
} // end double MatPat::dCdqlocal(

double MatPat::dCdsmlocmm(
        bool envt_deriv
        ,bool const envt_i
        ,bool const envt_j
        ,Sex const sex_j)
{
    double cval = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // information from a male who mates at the focal patch
        // but who is from a nonlocal patch himself
        double sum_male_nonlocal_mating = 0.0;

        for (int remote_envt_i = 0; remote_envt_i < 2; ++remote_envt_i)
        {
            sum_male_nonlocal_mating += 
                f(envt_i) * px(z_off_i, envt_i, remote_envt_i);

        }

        cval += // philopatric contribution to competing juveniles
            n[Female] * (1.0 - d[sex_j]) * (
                dFdsmlocmm(envt_deriv, envt_i) * 
                ( l * px(z_off_i, envt_i, envt_i)
                  + (1.0 - l) * sum_male_nonlocal_mating
                ) 
                +
                F(envt_i) *
                    l * dpxdsm(envt_deriv, z_off_i, envt_i, envt_i)
            ) * omega(z_off_i, envt_j);
    }

    return(cval);
} // end MatPat::dCdsmlocmm(

double MatPat::dCdsflocff(
        bool envt_deriv
        ,bool const envt_i
        ,bool const envt_j
        ,Sex const sex_j)
{
    double cval = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // information from a male who mates at the focal patch
        // but who is from a nonlocal patch himself
        double sum_male_nonlocal_mating = 0.0;

        double d_sum_male_nonlocal_mating_dsf = 0.0;

        for (int remote_envt_i = 0; remote_envt_i < 2; ++remote_envt_i)
        {
            sum_male_nonlocal_mating += 
                f(envt_i) * px(z_off_i, envt_i, remote_envt_i);

            d_sum_male_nonlocal_mating_dsf +=
                f(envt_i) * dpxdsf(envt_deriv, z_off_i, envt_i, remote_envt_i);
        }

        cval += // philopatric contribution to competing juveniles
            n[Female] * (1.0 - d[sex_j]) * (
                dFdsflocff(envt_deriv, envt_i) * 
                ( l * px(z_off_i, envt_i, envt_i)
                  + (1.0 - l) * sum_male_nonlocal_mating
                ) 
                +
                F(envt_i) *
                (
                    l * dpxdsf(envt_deriv, z_off_i, envt_i, envt_i)
                    + (1.0 - l) * d_sum_male_nonlocal_mating_dsf
                )
            ) * omega(z_off_i, envt_j);
    }

    return(cval);
}// end MatPat::dCdsflocff

// total number of competiting juveniles
// of sex_j
// in a patch that transitions
// from envt_i to envt_j,
double MatPat::C(
        bool const envt_i
        ,bool const envt_j
        ,Sex const sex_j)
{
    double cval = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // information from a male who mates at the focal patch
        // but who is from a nonlocal patch himself
        double sum_male_nonlocal_mating = 0.0;

        // information from a male who mates at a remote patch
        // but who is from another patch himself
        double sum_male_nonlocal_mating_at_remote_patch = 0.0;

        double number_immigrant_juveniles = 0.0;

        for (int remote_envt_i = 0; remote_envt_i < 2; ++remote_envt_i)
        {
            sum_male_nonlocal_mating += 
                f(envt_i) * px(z_off_i, envt_i, remote_envt_i);

            // calculate contribution of mal
            sum_male_nonlocal_mating_at_remote_patch = 0.0;

            for (int remote_envt_j = 0; remote_envt_j < 2; ++remote_envt_j)
            {
                sum_male_nonlocal_mating_at_remote_patch +=
                    f(remote_envt_j) * px(z_off_i, remote_envt_i, remote_envt_j);
            }

            number_immigrant_juveniles += f(remote_envt_i) *
                F(remote_envt_i) * 
                (l * px(z_off_i, remote_envt_i, remote_envt_i)
                 + (1.0 - l) * sum_male_nonlocal_mating_at_remote_patch);
        }

        cval += // philopatric contribution to competing juveniles
            n[Female] * (1.0 - d[sex_j]) * 
            F(envt_i) * 
            ( l * px(z_off_i, envt_i, envt_i)
              + (1.0 - l) * sum_male_nonlocal_mating
            ) * omega(z_off_i, envt_j)
            +
            n[Female] * d[sex_j] * number_immigrant_juveniles * omega(z_off_i, envt_j);
    }

    return(cval);
}

double MatPat::f(bool const envt)
{
    return(sigma[!envt] / (sigma[0] + sigma[1]));
}


// juvenile survival 
double MatPat::omega(
        bool const offspring_phenotype
        ,bool const envt_j)
{
    if (envt_j != offspring_phenotype)
    {
        return(1.0 - c[envt_j]);
    }

    return(1.0);
}


double MatPat::envt_switch(
        bool const envt_i
        ,bool const envt_j)
{
    if (envt_i == envt_j)
    {
        return(1.0 - sigma[envt_i]);
    }

    return(sigma[envt_i]);
}

// local fecundity
double MatPat::F(
        bool const envt_i
)
{
    // calculate the fecundity function which is 
    // 1/B, where B is the average cost
     double denominator = 0.0;

     for (size_t z_off = 0; z_off < 2; ++z_off)
     {
         // local mating
         denominator += y[z_off] * 
             l * px(z_off, envt_i, envt_i);

         for (size_t envt_remote = 0; envt_remote < 2; ++envt_remote)
         {
             denominator += y[z_off] * (1.0 - l)
                 * f(envt_remote) * px(z_off, envt_i, envt_remote);
         }
     }

     return(1.0 / (denominator * n[Female]));
}


double MatPat::dpxdq(
    int const q_deriv 
    ,bool const offspring_phenotype
    ,bool const envt_female
    ,bool const envt_male)
{
    assert(q_deriv >= 0);
    assert(q_deriv <= 2);

    double psignal_mom = (1.0 - error[Female]) * sf[envt_female]
        + error[Female] * sf[!envt_female];

    double psignal_dad = (1.0 - error[Male]) * sm[envt_male]
         + error[Male] * sm[!envt_male];

    assert(psignal_mom >= 0.0);
    assert(psignal_mom <= 1.0);
    assert(psignal_dad >= 0.0);
    assert(psignal_dad <= 1.0);

    double dq[3] = {
        (1.0 - psignal_mom) * (1.0 - psignal_dad),
        ((1.0 - psignal_mom) * psignal_dad +
        psignal_mom * (1.0 - psignal_dad)),
        psignal_mom * psignal_dad
    };

    // prob of producing offspring phenotype z2 is requested
    if (offspring_phenotype)
    {
        return(-dq[q_deriv]);
    }

    return(dq[q_deriv]);
} // end double MatPat::dpxdq(

double MatPat::dpxdsm(
    bool const envt_deriv
    ,bool const offspring_phenotype
    ,bool const envt_female
    ,bool const envt_male)
{
    double psignal_mom = (1.0 - error[Female]) * sf[envt_female]
        + error[Female] * sf[!envt_female];
    
    assert(psignal_mom >= 0.0);
    assert(psignal_mom <= 1.0);

    double dpsignal_dad_dsm = 0.0;

    if (envt_deriv == envt_male)
    {
        dpsignal_dad_dsm = (1.0 - error[Male]);
    }
    else
    {
        dpsignal_dad_dsm = error[Male];
    }
    
    double d_prob_offspring_z1_dsm = 
        psignal_mom * dpsignal_dad_dsm * q[2]
        + (psignal_mom * -dpsignal_dad_dsm +
        + (1.0 - psignal_mom) * dpsignal_dad_dsm) * q[1]
        + (1.0 - psignal_mom) * -dpsignal_dad_dsm * q[0];

    // prob of producing offspring phenotype z2 is requested
    if (offspring_phenotype)
    {
        return(-d_prob_offspring_z1_dsm);
    }

    return(d_prob_offspring_z1_dsm);
} // end MatPat::dpxdsm(


// derivative of the phenotype determination function px
// with respect to the female signal sf
double MatPat::dpxdsf(
    bool const envt_deriv
    ,bool const offspring_phenotype
    ,bool const envt_female
    ,bool const envt_male)
{
    double dpsignal_mom_dsf = 0.0;

    if (envt_deriv == envt_female)
    {
        dpsignal_mom_dsf = (1.0 - error[Female]);
    }
    else
    {
        dpsignal_mom_dsf = error[Female];
    }

    double psignal_dad = (1.0 - error[Male]) * sm[envt_male]
         + error[Male] * sm[!envt_male];

    assert(psignal_dad >= 0.0);
    assert(psignal_dad <= 1.0);

    double d_prob_offspring_z1_dsf = 
        dpsignal_mom_dsf * psignal_dad * q[2]
        + (dpsignal_mom_dsf * (1.0 - psignal_dad) +
        + -dpsignal_mom_dsf * psignal_dad) * q[1]
        + -dpsignal_mom_dsf * (1.0 - psignal_dad) * q[0];

    // prob of producing offspring phenotype z2 is requested
    if (offspring_phenotype)
    {
        return(-d_prob_offspring_z1_dsf);
    }

    return(d_prob_offspring_z1_dsf);
} // end MatPat::dpxdsf

double MatPat::dFdsmlocmm(
        bool const envt_deriv // derivative of sf(envt_deriv)
        ,bool const envt_i
        )
{
    double val = 0.0;

    double denom = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            denom += n[Female] * y[z_off_i] * (1.0 - l) * f(envt_remote) *
                px(z_off_i, envt_i, envt_remote);
        }

        val += n[Female] * y[z_off_i] * l * 
            dpxdsm(envt_deriv, z_off_i, envt_i, envt_i);

        denom += n[Female] * y[z_off_i] * l * px(z_off_i, envt_i, envt_i);
    }

    return(-val/(denom * denom));
} // end MatPat::dFdsmlocmm

double MatPat::dFdsflocff(
        bool const envt_deriv // derivative of sf(envt_deriv)
        ,bool const envt_i
        )
{
    double val = 0.0;

    double denom = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            val += n[Female] * y[z_off_i] * (1.0 - l) * f(envt_remote) * 
                dpxdsf(envt_deriv, z_off_i, envt_i, envt_remote);

            denom += n[Female] * y[z_off_i] * (1.0 - l) * f(envt_remote) *
                px(z_off_i, envt_i, envt_remote);
        }

        val += n[Female] * y[z_off_i] * l * 
            dpxdsf(envt_deriv, z_off_i, envt_i, envt_i);

        denom += n[Female] * y[z_off_i] * l * px(z_off_i, envt_i, envt_i);
    }

    return(-val/(denom * denom));
}

// derivative of F for a male who sires from a remote patch
double MatPat::dFdqremote(
        int const combn_i_deriv // derivative of sf(envt_deriv)
        ,bool const envt_i
        )
{
    double val = 0.0;

    double denom = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            val += n[Female] * y[z_off_i] * (1.0 - l) * f(envt_remote) * 
                dpxdq(combn_i_deriv, z_off_i, envt_i, envt_remote);

            denom += n[Female] * y[z_off_i] * (1.0 - l) * f(envt_remote) *
                px(z_off_i, envt_i, envt_remote);
        }

        denom += n[Female] * y[z_off_i] * l * px(z_off_i, envt_i, envt_i);
    }

    return(-val/(denom * denom));
} // double MatPat::dFdqremote(

double MatPat::dFdqloc(
        int const combn_i_deriv // derivative of sf(envt_deriv)
        ,bool const envt_i
        )
{
    double val = 0.0;

    double denom = 0.0;

    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            denom += n[Female] * y[z_off_i] * (1.0 - l) * f(envt_remote) *
                px(z_off_i, envt_i, envt_remote);
        }

        val += n[Female] * y[z_off_i] * l * 
            dpxdq(combn_i_deriv, z_off_i, envt_i, envt_i);

        denom += n[Female] * y[z_off_i] * l * px(z_off_i, envt_i, envt_i);
    }

    return(-val/(denom * denom));
} // double MatPat::dFdqloc(

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt sflocmf[envt_deriv]
double MatPat::dWdsflocmf(
        int dim_i
        ,int dim_j
        ,bool envt_deriv
        )// which of the two envt'al signals to take derivs over
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;

    Sex sex_current = dim_j > 1 ? Male : Female;


    if (sex_current == Female)
    {
        return(0.0);
    }

    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        double locally_sired_remote_established = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            locally_sired_remote_established += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);
        }

        // w = (1/2) * nf/nm * l * F' * (
        //              px * omega * (local + remote)
        //
        //  + (1/2) * nf/nm * l * F  (
        //              px' * omega * (local + remote)
        //              )
        //  + (1/2) * nf/nm * l * F  (
        //              px * omega * (local' + remote')
        //              )
        val += 0.5 * n[Female] / n[Male] * l * 
                    (
                         dFdsflocff(envt_deriv, envt_current) *
                                        px(z_off_i, envt_current, envt_current) 
                        + F(envt_current) * 
                                        dpxdsf(envt_deriv, z_off_i, envt_current, envt_current) 
                    )
            * omega(z_off_i, envt_future) *
                    ( 
                        n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future) 
                            * envt_switch(envt_current, envt_future)
                        + locally_sired_remote_established
                    );

        val += 0.5 * n[Female] / n[Male] * l * F(envt_current)
                    * px(z_off_i, envt_current, envt_current)
                    * omega(z_off_i, envt_future)
                    * (
                            n[sex_future] * (1.0 - d[sex_future]) * 
                            envt_switch(envt_current,envt_future) *
                                -dCdsflocff(envt_deriv, envt_current, envt_future, sex_future) /
                                pow(C(envt_current, envt_future, sex_future),2)
                    );
    }
    return(val);
}

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt dWdsmlocfm[envt_deriv]
double MatPat::dWdsmlocfm(
        int dim_i // future envt
        ,int dim_j // current envt
        ,bool envt_deriv // which of the two envt'al signals to take derivs over
    ) 
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;

    if (sex_current == Male)
    {
        return(0.0);
    }

    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        double sum_remote_siring = 0.0;
        double sum_establishment_remote = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_remote_siring += f(envt_remote) * 
                px(z_off_i, envt_current, envt_remote);

            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);
        }

        // Let w = F(smlocfm) k(smlocfm) E(smlocfm), where F is fecundity and k are 
        // probabilities of local vs nonlocal mating
        // and E is the probability of establishment
        //
        // then dw/dsmlocfm = 
        // F'(sm) k(sm) E(sm) +
        // F(sm) k'(sm) E(sm) +
        // F(sm) k(sm) E'(sm) +
        val += 0.5 * dFdsmlocmm(envt_deriv, envt_current) *
                        (
                            l * px(z_off_i, envt_current, envt_current)
                            + (1.0 - l) * sum_remote_siring
                        )
                   * omega(z_off_i, envt_future) * 
                        (
                            envt_switch(envt_current, envt_future) *
                            n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote
                        );

        val += 0.5 * F(envt_current) *
                        l * dpxdsm(envt_deriv, z_off_i, envt_current, envt_current) *
                    omega(z_off_i, envt_future) * 
                        (
                            envt_switch(envt_current, envt_future) *
                            n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote
                        );


        val += 0.5 * F(envt_current) *
                        (
                            l * px(z_off_i, envt_current, envt_current)
                            + (1.0 - l) * sum_remote_siring
                        )
                * omega(z_off_i, envt_future) *
                        (
                            envt_switch(envt_current, envt_future) 
                                * n[sex_future] * (1.0 - d[sex_future]) * 
                                    -dCdsmlocmm(
                                        envt_deriv,
                                        envt_current,
                                        envt_future,
                                        sex_future) /
                            pow(C(envt_current, envt_future, sex_future),2)
                        );
    }
    
    return(val);
}

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt sflocff[envt_deriv]
double MatPat::dWdsflocff(
        int dim_i // future envt
        ,int dim_j // current envt
        ,bool envt_deriv // which of the two envt'al signals to take derivs over
    ) 
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;

    if (sex_current == Male)
    {
        return(0.0);
    }

    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        double sum_remote_siring = 0.0;
        double sum_establishment_remote = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_remote_siring += f(envt_remote) * 
                px(z_off_i, envt_current, envt_remote);

            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);
        }

        // Let w = F(sflocff) k(sffoc) E(sflocff), where F is fecundity and k are 
        // probabilities of local vs nonlocal mating
        // and E is the probability of establishment
        //
        // then dw/dsflocff = F'(sf) k(sf) E(sf) + F(sf) k(sf) E'(sflocff)
        val += 0.5 * dFdsflocff(envt_deriv, envt_current) *
                        (
                            l * px(z_off_i, envt_current, envt_current)
                            + (1.0 - l) * sum_remote_siring
                        )
                   * omega(z_off_i, envt_future) * 
                        (
                            envt_switch(envt_current, envt_future) *
                            n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote
                        );


        val += 0.5 * F(envt_current) *
                        (
                            l * px(z_off_i, envt_current, envt_current)
                            + (1.0 - l) * sum_remote_siring
                        )
                * omega(z_off_i, envt_future) *
                        (
                            envt_switch(envt_current, envt_future) 
                                * n[sex_future] * (1.0 - d[sex_future]) * 
                                    -dCdsflocff(
                                        envt_deriv,
                                        envt_current,
                                        envt_future,
                                        sex_future) /
                            pow(C(envt_current, envt_future, sex_future),2)
                        );
    }
    
    return(val);
}


// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt smlocmm[envt_i]
double MatPat::dWdsmlocmm(
        int dim_i // future
        ,int dim_j // current
        ,bool envt_deriv // which of the two envt'al signals to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;
   
    // if female, she cannot express this trait, hence 0.0
    if (sex_current == Female)
    {
        return(0.0);
    }

    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double sum_establishment_remote = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);
        }

        // Let w = F(smlocfm) k(smfoc) E(smlocfm), where F is fecundity and k are 
        // probabilities of local vs nonlocal mating
        // and E is the probability of establishment
        //
        // then dw/dsmfoc = F(sm) k'(sm) E(sm)
        //
        //
        val += 0.5 * n[Female]/n[Male] * l * dFdsmlocmm(envt_deriv,envt_current) 
                * px(z_off_i, envt_current, envt_current)
                * omega(z_off_i, envt_future)
            * (
                    envt_switch(envt_current, envt_future) * 
                    n[sex_future] * (1.0 - d[sex_future]) /
                    C(envt_current, envt_future, sex_future)
                    +
                    sum_establishment_remote
            );
            //+ 0.5 * n[Female]/n[Male] * (1.0 - l) * sum_nonlocal_sirings;
            //
        val += 0.5 * n[Female]/n[Male] * l * F(envt_current)
                * px(z_off_i, envt_current, envt_current)
                * omega(z_off_i, envt_future)
                *    envt_switch(envt_current, envt_future) * 
                    n[sex_future] * (1.0 - d[sex_future]) *
                        -dCdsmlocmm(
                            envt_deriv
                            ,envt_current
                            ,envt_future
                            ,sex_future) /
                    pow(C(envt_current, envt_future, sex_future),2);
    }
    
    return(val);
}

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt sffoc[envt_i]
double MatPat::dWdsffoc(
        int dim_i // future
        ,int dim_j // current
        ,bool envt_deriv // which of the two envt'al signals to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;
    
    if (sex_current == Male)
    {
        return(0.0);
    }

    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double sum_remote_siring = 0.0;
        double d_sum_remote_siring_dsffoc = 0.0;
        double sum_establishment_remote = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_remote_siring += f(envt_remote) * 
                px(z_off_i, envt_current, envt_remote);

            d_sum_remote_siring_dsffoc += f(envt_remote) * 
                dpxdsf(envt_deriv, z_off_i, envt_current, envt_remote);

            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);
        }

        // Let w = F(sflocff) k(sffoc) E(sflocff), where F is fecundity and k are 
        // probabilities of local vs nonlocal mating
        // and E is the probability of establishment
        //
        // then dw/dsffoc = F(sf) k'(sf) E(sf)
        //
        //
        val += 0.5 * F(envt_current) *
                        (
                            l * dpxdsf(envt_deriv, z_off_i, envt_current, envt_current)
                            + (1.0 - l) * d_sum_remote_siring_dsffoc
                        )
               * omega(z_off_i, envt_future) * 
                        (
                            envt_switch(envt_current, envt_future) *
                                n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote
                        );
    }
    
    return(val);
}


// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt smfoc[envt_i]
double MatPat::dWdsmfoc(
        int dim_i // future
        ,int dim_j // current
        ,bool envt_deriv // which of the two envt'al signals to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;
   
    // if female, she cannot express this trait, hence 0.0
    if (sex_current == Female)
    {
        return(0.0);
    }

    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double sum_establishment_remote = 0.0;
        double d_sum_nonlocal_sirings_dsm = 0.0;
        double sum_nonlocal_sirings_establish_remote = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_nonlocal_sirings_establish_remote = 0.0;

            for (int envt_remote_2 = 0; envt_remote_2 < 2; ++envt_remote_2)
            {
                sum_nonlocal_sirings_establish_remote += 
                    f(envt_remote_2) * envt_switch(envt_remote_2, envt_future)
                    * n[sex_future] * d[sex_future] / 
                    C(envt_remote_2, envt_future, sex_future);
            }


            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);

            d_sum_nonlocal_sirings_dsm += f(envt_remote)
                * F(envt_remote) * dpxdsm(envt_deriv, z_off_i, envt_remote, envt_current)
                * omega(z_off_i, envt_future)
                * (
                        envt_switch(envt_remote, envt_future) * n[sex_future] * 
                        (1.0 - d[sex_future]) / C(envt_remote, envt_future, sex_future)
                        +
                        sum_nonlocal_sirings_establish_remote
                );
        }

        // Let w = F(smlocfm) k(smfoc) E(smlocfm), where F is fecundity and k are 
        // probabilities of local vs nonlocal mating
        // and E is the probability of establishment
        //
        // then dw/dsmfoc = F(sm) k'(sm) E(sm)
        //
        //
        val += 0.5 * n[Female]/n[Male] * l * F(envt_current) 
            * dpxdsm(envt_deriv, z_off_i, envt_current, envt_current)
            * omega(z_off_i, envt_future)
            * (
                    envt_switch(envt_current, envt_future) * 
                    n[sex_future] * (1.0 - d[sex_future]) /
                    C(envt_current, envt_future, sex_future)
                    +
                    sum_establishment_remote
            )
            + 0.5 * n[Female]/n[Male] * (1.0 - l) * d_sum_nonlocal_sirings_dsm;
    }
    
    return(val);
} // end double MatPat::dWdsmfoc(


// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt qlocf
double MatPat::dWdqlocf(
        int dim_i // future
        ,int dim_j // current
        ,int combn_i_deriv // which of the 4 responsiveness traits to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;

    // this only pertains to fitness accrued through Females
    if (sex_current == Male)
    {
        return(0.0);
    }
   
    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double sum_establishment_remote = 0.0;
        double sum_nonlocal_mating = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);

            sum_nonlocal_mating += f(envt_remote) * px(z_off_i, envt_current, envt_remote);
        }

        val += 0.5 * dFdqloc(combn_i_deriv, envt_current) 
                    * (
                            l * px(z_off_i, envt_current, envt_current)
                            + (1.0 - l) * sum_nonlocal_mating
                    )
                    * omega(z_off_i, envt_future)
                    * (
                            envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote
                    );

        val += 0.5 * F(envt_current)
                    * (
                            l * px(z_off_i, envt_current, envt_current)
                            + (1.0 - l) * sum_nonlocal_mating
                    )
                    * omega(z_off_i, envt_future)
                    * envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) *
                            -dCdqlocal(combn_i_deriv, envt_current, envt_future, sex_future) /
                            pow(C(envt_current, envt_future, sex_future),2);
    }
    
    return(val);
} // end double MatPat::dWdqlocf

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt qlocfremote
double MatPat::dWdqlocfremote(
        int dim_i // future
        ,int dim_j // current
        ,int combn_i_deriv // which of the 4 responsiveness traits to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;

    // this only pertains to fitness accrued through Females
    if (sex_current == Male)
    {
        return(0.0);
    }
   
    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double sum_establishment_remote = 0.0;
        double sum_nonlocal_mating = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);

            sum_nonlocal_mating += f(envt_remote) * px(z_off_i, envt_current, envt_remote);
        }

        val += 0.5 * dFdqremote(combn_i_deriv, envt_current) 
                    * (
                            l * px(z_off_i, envt_current, envt_current)
                            + (1.0 - l) * sum_nonlocal_mating
                    )
                    * omega(z_off_i, envt_future)
                    * (
                            envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote
                    );

        val += 0.5 * F(envt_current)
                    * (
                            l * px(z_off_i, envt_current, envt_current)
                            + (1.0 - l) * sum_nonlocal_mating
                    )
                    * omega(z_off_i, envt_future)
                    *   envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) *
                            -dCdqremote(combn_i_deriv, envt_current, envt_future, sex_future) /
                            pow(C(envt_current, envt_future, sex_future),2);
    }
    
    return(val);
} // end double MatPat::dWdqlocfremote

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt qfocfremote
double MatPat::dWdqfocfremote(
        int dim_i // future
        ,int dim_j // current
        ,int combn_i_deriv // which of the 4 responsiveness traits to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;

    // this only pertains to fitness accrued through Females
    if (sex_current == Male)
    {
        return(0.0);
    }
   
    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double sum_establishment_remote = 0.0;
        double d_sum_nonlocal_mating_qfocfremote = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);

            d_sum_nonlocal_mating_qfocfremote += 
                f(envt_remote) * dpxdq(combn_i_deriv, z_off_i, envt_current, envt_remote);
        }

        val += 0.5 * F(envt_current) 
                    * (1.0 - l) * d_sum_nonlocal_mating_qfocfremote
                    * omega(z_off_i, envt_future)
                    * (
                            envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote
                    );
    }

    return(val);
} // end double MatPat::dWdqfocfremote

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt qfocm
double MatPat::dWdqfocmremote(
        int dim_i // future
        ,int dim_j // current
        ,int combn_i_deriv // which of the 4 responsiveness traits to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;

    // this only pertains to fitness accrued through Males
    if (sex_current == Female)
    {
        return(0.0);
    }
   
    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double d_sum_sirings_remote = 0.0;

        double sum_sirings_remote_dispersed = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_sirings_remote_dispersed = 0.0;

            for (int envt_remote_2 = 0; envt_remote_2 < 2; ++envt_remote_2)
            {
                sum_sirings_remote_dispersed += f(envt_remote_2)
                    * envt_switch(envt_remote_2, envt_future)
                    * n[sex_future] * d[sex_future] /
                        C(envt_remote_2, envt_future, sex_future);
            }

            d_sum_sirings_remote += f(envt_remote)
                * F(envt_remote)
                * dpxdq(combn_i_deriv, z_off_i, envt_remote, envt_current)
                * omega(z_off_i, envt_future)
                * (
                        envt_switch(envt_remote, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) /
                        C(envt_remote, envt_future, sex_future)
                        +
                        sum_sirings_remote_dispersed
                );
        }

        val += 0.5 * (1.0 - l) * n[Female]/n[Male] * d_sum_sirings_remote;
    }
    
    return(val);
} // end double MatPat::dWdqfocmremote(

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt qfocm
double MatPat::dWdqlocmremote(
        int dim_i // future
        ,int dim_j // current
        ,int combn_i_deriv // which of the 4 responsiveness traits to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;

    // this only pertains to fitness accrued through Males
    if (sex_current == Female)
    {
        return(0.0);
    }
   
    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double sum_establishment_remote = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);
        }

        val += 0.5 * l * n[Female]/n[Male] * dFdqremote(combn_i_deriv, envt_current) 
                    * px(z_off_i, envt_current, envt_current)
                    * omega(z_off_i, envt_future)
                    * (
                            envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote
                    );

        val += 0.5 * l * n[Female]/n[Male] * F(envt_current)
                    * px(z_off_i, envt_current, envt_current)
                    * omega(z_off_i, envt_future)
                    * envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) *
                            -dCdqremote(combn_i_deriv, envt_current, envt_future, sex_future) /
                            pow(C(envt_current, envt_future, sex_future),2);
    }
    
    return(val);
} // end double MatPat::dWdqlocmremote(

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt qfocm
double MatPat::dWdqlocm(
        int dim_i // future
        ,int dim_j // current
        ,int combn_i_deriv // which of the 4 responsiveness traits to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;

    // this only pertains to fitness accrued through Males
    if (sex_current == Female)
    {
        return(0.0);
    }
   
    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double sum_establishment_remote = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);
        }

        val += 0.5 * l * n[Female]/n[Male] * dFdqloc(combn_i_deriv, envt_current) 
                    * px(z_off_i, envt_current, envt_current)
                    * omega(z_off_i, envt_future)
                    * (
                            envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote
                    );

        val += 0.5 * l * n[Female]/n[Male] * F(envt_current)
                    * px(z_off_i, envt_current, envt_current)
                    * omega(z_off_i, envt_future)
                    * envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) *
                            -dCdqlocal(combn_i_deriv, envt_current, envt_future, sex_future) /
                            pow(C(envt_current, envt_future, sex_future),2);
    }
    
    return(val);
} // end double MatPat::dWdqlocm(

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt qfocm
double MatPat::dWdqfocm(
        int dim_i // future
        ,int dim_j // current
        ,int combn_i_deriv // which of the 4 responsiveness traits to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;

    // this only pertains to fitness accrued through Males
    if (sex_current == Female)
    {
        return(0.0);
    }
   
    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double sum_establishment_remote = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_establishment_remote += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);
        }

        val += 0.5 * n[Female]/n[Male] * l * F(envt_current) 
                    * dpxdq(combn_i_deriv, z_off_i, envt_current, envt_current)
                    * omega(z_off_i, envt_future)
                    * (
                            envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote
                    );
    }
    
    return(val);
} // end double MatPat::dWdqfocm(

// derivative of the transition matrix B
// element b_(dim_i, dim_j) wrt qfocf
double MatPat::dWdqfocf(
        int dim_i // future
        ,int dim_j // current
        ,int combn_i_deriv // which of the 4 responsiveness traits to take derivs over
        )
{
    double val = 0.0;

    // infer future environment 
    // from matrix dimensions
    // 0, 2: envt 0
    // 1, 3: envt 1
    int envt_future = dim_i % 2;

    // infer current environment
    // from matrix dimensions
    int envt_current = dim_j % 2;

    // infer future sex from matrix dimensions
    // 0, 1: female
    // 2, 3: male
    Sex sex_future = dim_i > 1 ? Male : Female;
    
    Sex sex_current = dim_j > 1 ? Male : Female;

    // this only pertains to fitness accrued through Females
    if (sex_current == Male)
    {
        return(0.0);
    }
   
    // calculate fitness by calculating the expected probability
    // that offspring of phenotype z_off_i establish themselves
    // locally or remotely
    for (int z_off_i = 0; z_off_i < 2; ++z_off_i)
    {
        // sum all fitness components that 
        // involve remote environments
        double sum_establishment_remote_f = 0.0;
        double nonlocal_px = 0.0;

        for (int envt_remote = 0; envt_remote < 2; ++envt_remote)
        {
            sum_establishment_remote_f += f(envt_remote)
                * envt_switch(envt_remote, envt_future) * n[sex_future] * d[sex_future] /
                        C(envt_remote, envt_future, sex_future);

            nonlocal_px += f(envt_remote) * px(z_off_i, envt_current, envt_remote);
        }

        val += 0.5 * F(envt_current)
                    * l * dpxdq(combn_i_deriv, z_off_i,  envt_current, envt_current)
                    * omega(z_off_i, envt_future)
                    * (
                            envt_switch(envt_current, envt_future) * 
                            n[sex_future] * (1.0 - d[sex_future]) /
                            C(envt_current, envt_future, sex_future)
                            +
                            sum_establishment_remote_f
                    );
    }
    
    return(val);
} // end double MatPat::dWdqfocf

// selection gradient on female signal
double MatPat::dWdsf(bool const envt_i)
{
    double dwdsf = 0.0;

    // relatedness coefficients
    double rMfocf = 0.0; // relatedness between focal mother and herself
    double rMlocff = 0.0; // relatedness between focal mother and any local female breeder
    double rMlocmf = 0.0; // relatedness between focal mother and any local male breeder


    // aux var to get the maternal environment 
    // from the matrix dimensions
    bool envt_from_dim;

    // iterate through future
    // environment x sex combinations
    for (size_t dim_i = 0; dim_i < 4; ++dim_i)
    {
        // iterate through previous 
        // environment x sex combinations
        for (size_t dim_j = 0; dim_j < 4; ++dim_j)
        {
            // get the environment from the matrix dimensions
            envt_from_dim = dim_j % 2;

            rMfocf = 1.0;

            rMlocff = 1.0 / n[Female] + 
                (n[Female] - 1.0)/n[Female] * Q[Female][Female][envt_from_dim];

            rMlocmf = Q[Female][Male][envt_from_dim];

            dwdsf += v[dim_i] * u[dim_j] * 
                (
                    dWdsffoc(dim_i, dim_j, envt_i) * rMfocf
                    +
                    dWdsflocff(dim_i, dim_j, envt_i) * rMlocff
                    +
                    dWdsflocmf(dim_i, dim_j, envt_i) * rMlocmf
                );
        } // end for dim_j
    }

    return(dwdsf);
} // end double MatPat::dWdsf(bool const envt_i)


// selection gradient on male signal
double MatPat::dWdsm(bool const envt_i)
{
    double dwdsm = 0.0;

    // relatedness coefficients
    double rPfocm = 0.0; // relatedness between focal father and herself
    double rPlocmm = 0.0; // relatedness between focal father and any local male breeder
    double rPlocfm = 0.0; // relatedness between focal father and any local female breeder


    // aux var to get the maternal environment 
    // from the matrix dimensions
    bool envt_from_dim;

    // iterate through future
    // environment x sex combinations
    for (size_t dim_i = 0; dim_i < 4; ++dim_i)
    {
        // iterate through previous 
        // environment x sex combinations
        for (size_t dim_j = 0; dim_j < 4; ++dim_j)
        {
            // get the environment from the matrix dimensions
            envt_from_dim = dim_j % 2;

            rPfocm = 1.0;

            rPlocmm = 1.0 / n[Male] + 
                (n[Male] - 1.0)/n[Male] * Q[Male][Male][envt_from_dim];

            rPlocfm = Q[Male][Female][envt_from_dim];

            dwdsm += v[dim_i] * u[dim_j] * 
                (
                    dWdsmfoc(dim_i, dim_j, envt_i) * rPfocm
                    +
                    dWdsmlocmm(dim_i, dim_j, envt_i) * rPlocmm
                    +
                    dWdsmlocfm(dim_i, dim_j, envt_i) * rPlocfm
                );
        } // end for dim_j
    }

    return(dwdsm);
}


// selection gradient on offspring responsiveness
double MatPat::dWdq(int const combn_i)
{
    double dwdq = 0.0;

    // relatedness coefficients
    double rfocf = 0.0; // relatedness between focal offspring and itself
    double rfocf_remote = 0.0; // relatedness between focal offspring and any remotely sired offspring
    double rlocf = 0.0; // relatedness between focal father and any local female breeder
    double rlocf_remote = 0.0; // relatedness between focal father and any local female breeder

    double rfocm = 0.0; // relatedness between focal father and herself
    double rfocm_remote = 0.0; // relatedness between focal father and any local male breeder
    double rlocm = 0.0; // relatedness between focal father and any local female breeder
    double rlocm_remote = 0.0; // relatedness between focal father and any local female breeder

    // aux var to get the maternal environment 
    // from the matrix dimensions
    bool envt_from_dim;

    // iterate through future
    // environment x sex combinations
    for (size_t dim_i = 0; dim_i < 4; ++dim_i)
    {
        // iterate through previous 
        // environment x sex combinations
        for (size_t dim_j = 0; dim_j < 4; ++dim_j)
        {
            // get the environment from the matrix dimensions
            envt_from_dim = dim_j % 2;

            rfocf = 1.0;

            rfocf_remote = 1.0;

            rlocf = 0.5 * (1.0/n[Female] + (n[Female] - 1.0) / n[Female] * 
                    Q[Female][Female][envt_from_dim])
                    + 0.5 * Q[Female][Male][envt_from_dim];

            rlocf_remote = 0.5 * (1.0/n[Female] + 
                    (n[Female] - 1.0)/n[Female] * Q[Female][Female][envt_from_dim]);


            rfocm = 1.0;

            rfocm_remote = 1.0;

            rlocm = 0.5 * (1.0/n[Male] + (n[Male] - 1.0) /n[Male] * Q[Male][Male][envt_from_dim])
                + 0.5 * Q[Female][Male][envt_from_dim];

            rlocm_remote = 0.5 * Q[Female][Male][envt_from_dim];

            dwdq += v[dim_i] * u[dim_j] * 
                (
                     dWdqfocf(dim_i, dim_j, combn_i) * rfocf
                    +
                  dWdqfocfremote(dim_i, dim_j, combn_i) * rfocf_remote
                  +
                   dWdqlocf(dim_i, dim_j, combn_i) * rlocf
                  +
                  dWdqlocfremote(dim_i, dim_j, combn_i) * rlocf_remote
                  +
                   dWdqfocm(dim_i, dim_j, combn_i) * rfocm
                  +
                  dWdqfocmremote(dim_i, dim_j, combn_i) * rfocm_remote
                  +
                  dWdqlocm(dim_i, dim_j, combn_i) * rlocm
                    +
                    dWdqlocmremote(dim_i, dim_j, combn_i) * rlocm_remote
                );

        } // end for dim_j
    }

    return(dwdq);
}
