#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "bramauxiliary.h"


using namespace std;

string filename("iter_dispersal_status");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

struct rparams
{
    // VARS
};

double bound01(double val)
{
    val = val < 0 ? 0.0001 : val > 1.0 ? 0.9999 : val;

    return(val);
}

double bound0(double val)
{
    val = val < 0 ? 0.0001 : val;

    return(val);
}


// recursions of all the patch frequencies
int psys_recur(void *params, gsl_vector *f)
{
    //VARFUNC
    
    
    // PATCHRECUR
    
    return GSL_SUCCESS;
}


void reproductive_values(void *params, gsl_vector *f)
{
    // VARFUNC


    // REPVALS
}

void selgrads(void *params, gsl_vector *f)
{
    // VARFUNC

    // SELGRADS

}


void write_params(void *params)
{
    // VARFUNC


    // WRITEPARAMS
}


void write_data(void *params, int time)
{
    // VARFUNC


    if (time < 0)
    {
        // HEADERWRT
    }

    // WRITEDATA
}


int main (int argc, char **argv)
{
    int max_iter = atoi(argv[1]);

    // initialize the vectors that contain the variables
    // functions solve for
    //
    // vector for patch frequencies
    gsl_vector *x_p = gsl_vector_alloc(8);
    
    // vector for reproductive values
    gsl_vector *x_v = gsl_vector_alloc(8);

    // vector for selection gradients
    gsl_vector *x_selgrad = gsl_vector_alloc(4);

    // initialize the struct with parameters
    struct rparams paramstruct; 
  
    // initialize command line argument things
    // see generate_cpp.py 
    // ARGINIT
   

    // ranges for cycling: sometimes numerical iterations
    // won't resolve as solutions slightly cycle around
    // the final value. To see whether cyling behaviour
    // is repetitive over time, however, we create these
    // vectors that stores a series of past values of 
    // the values for the switching rates and compares
    // those to values that are currently found
    gsl_vector *p1imm_range = gsl_vector_alloc(10);
    gsl_vector *p2imm_range = gsl_vector_alloc(10);
    gsl_vector *p1phil_range = gsl_vector_alloc(10);
    gsl_vector *p2phil_range = gsl_vector_alloc(10);

    for (int ik = 0; ik < 10; ++ik)
    {
        // initialize the vector
        gsl_vector_set(p1imm_range, ik, 0);        
        gsl_vector_set(p2imm_range, ik, 0);        
        gsl_vector_set(p1phil_range, ik, 0);        
        gsl_vector_set(p2phil_range, ik, 0);        
    }

    // write the initial data set
    write_data(&paramstruct,-1);

    // iterate
    int iter;
    for (iter = 0; iter < max_iter ; ++iter)
    {
        // patch frequencies
        // cout << "patch" << endl;
        psys_recur(&paramstruct, x_p);

        paramstruct.f1aimm = gsl_vector_get(x_p,0);
        paramstruct.f2aimm = gsl_vector_get(x_p,1);
        paramstruct.f1mimm = gsl_vector_get(x_p,2);
        paramstruct.f2mimm = gsl_vector_get(x_p,3);
        paramstruct.f1aphil = gsl_vector_get(x_p,4);
        paramstruct.f2aphil = gsl_vector_get(x_p,5);
        paramstruct.f1mphil = gsl_vector_get(x_p,6);
        paramstruct.f2mphil = gsl_vector_get(x_p,7);

        assert(isnan(paramstruct.f1aimm) == 0);
        assert(isnan(paramstruct.f2aimm) == 0);
        assert(isnan(paramstruct.f1mimm) == 0);
        assert(isnan(paramstruct.f2mimm) == 0);
        assert(isnan(paramstruct.f1aphil) == 0);
        assert(isnan(paramstruct.f2aphil) == 0);
        assert(isnan(paramstruct.f1mphil) == 0);
        assert(isnan(paramstruct.f2mphil) == 0);

        // reproductive values
        // cout << "rval" << endl;
        reproductive_values(&paramstruct, x_v);
        
        paramstruct.v1aimm = gsl_vector_get(x_v, 0);
        paramstruct.v2aimm = gsl_vector_get(x_v, 1);
        paramstruct.v1mimm = gsl_vector_get(x_v, 2);
        paramstruct.v2mimm = gsl_vector_get(x_v, 3);
        paramstruct.v1aphil = gsl_vector_get(x_v, 4);
        paramstruct.v2aphil = gsl_vector_get(x_v, 5);
        paramstruct.v1mphil = gsl_vector_get(x_v, 6);
        paramstruct.v2mphil = gsl_vector_get(x_v, 7);
        
        assert(isnan(paramstruct.v1aimm) == 0);
        assert(isnan(paramstruct.v2aimm) == 0);
        assert(isnan(paramstruct.v1mimm) == 0);
        assert(isnan(paramstruct.v2mimm) == 0);
        assert(isnan(paramstruct.v1aphil) == 0);
        assert(isnan(paramstruct.v2aphil) == 0);
        assert(isnan(paramstruct.v1mphil) == 0);
        assert(isnan(paramstruct.v2mphil) == 0);

        // selection gradients
        // cout << "selgrad" << endl;
        selgrads(&paramstruct, x_selgrad);

        bool condition_p1imm = fabs(paramstruct.p1imm - gsl_vector_get(x_selgrad, 0)) < 1e-10; 
        bool condition_p2imm = fabs(paramstruct.p2imm - gsl_vector_get(x_selgrad, 0)) < 1e-10; 
        bool condition_p1phil = (fabs(paramstruct.p1phil - gsl_vector_get(x_selgrad, 0)) < 1e-10 || paramstruct.p1phil >= 0.999) || paramstruct.p1phil <= 0.001;
        bool condition_p2phil = (fabs(paramstruct.p2phil - gsl_vector_get(x_selgrad, 0)) < 1e-10 || paramstruct.p2phil >= 0.999) || paramstruct.p2phil <= 0.001;

        if (condition_p1imm && condition_p2imm && condition_p1phil && condition_p2phil)
        {
            paramstruct.p1imm = gsl_vector_get(x_selgrad, 0);
            paramstruct.p1phil = gsl_vector_get(x_selgrad, 1);
            paramstruct.p2imm = gsl_vector_get(x_selgrad, 2);
            paramstruct.p2phil = gsl_vector_get(x_selgrad, 3);

            break;
        }
        paramstruct.p1imm = gsl_vector_get(x_selgrad, 0);
        paramstruct.p1phil = gsl_vector_get(x_selgrad, 1);
        paramstruct.p2imm = gsl_vector_get(x_selgrad, 2);
        paramstruct.p2phil = gsl_vector_get(x_selgrad, 3);
        
        assert(isnan(paramstruct.p1imm) == 0);
        assert(isnan(paramstruct.p1phil) == 0);
        assert(isnan(paramstruct.p2imm) == 0);
        assert(isnan(paramstruct.p2phil) == 0);
        
        if (iter > 50000)
        {
            bool found_in_range1 = false;
            bool found_in_range2 = false;
            bool found_in_range3 = false;
            bool found_in_range4 = false;

            bool done = false;

            for (int ik = 0; ik < 10; ++ik)
            {
                if (!found_in_range1 && fabs(gsl_vector_get(p1imm_range, ik)-paramstruct.p1imm) < 1e-10)
                {
                    found_in_range1 = true;
                }

                if (!found_in_range2 && fabs(gsl_vector_get(p2imm_range, ik)-paramstruct.p2imm) < 1e-10)
                {
                    found_in_range2 = true;
                }
                
                if (!found_in_range3 && fabs(gsl_vector_get(p1phil_range, ik)-paramstruct.p1phil) < 1e-10)
                {
                    found_in_range3 = true;
                }
                if (!found_in_range4 && fabs(gsl_vector_get(p2phil_range, ik)-paramstruct.p2phil) < 1e-10)
                {
                    found_in_range4 = true;
                }
            }

            if (done)
            {
                break;
            }
        }


        for (int ik = 9; ik > 0; --ik)
        {
            gsl_vector_set(p1imm_range, ik, gsl_vector_get(p1imm_range, ik - 1));
            gsl_vector_set(p2imm_range, ik, gsl_vector_get(p2imm_range, ik - 1));
            gsl_vector_set(p1phil_range, ik, gsl_vector_get(p1phil_range, ik - 1));
            gsl_vector_set(p2phil_range, ik, gsl_vector_get(p2phil_range, ik - 1));
        }

        gsl_vector_set(p1imm_range, 0, paramstruct.p1imm);
        gsl_vector_set(p2imm_range, 0, paramstruct.p2imm);
        gsl_vector_set(p1phil_range, 0, paramstruct.p1phil);
        gsl_vector_set(p2phil_range, 0, paramstruct.p2phil);

        if (iter % 100 == 0)
        {
            write_data(&paramstruct,iter);
        }
    }

    write_data(&paramstruct,iter);
    write_params(&paramstruct);

    gsl_vector_free(x_p);
    gsl_vector_free(x_v);
    gsl_vector_free(x_selgrad);
}
