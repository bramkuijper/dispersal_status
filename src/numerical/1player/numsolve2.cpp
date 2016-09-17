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
    // 
double mum;
double mua;
double d;
double s1;
double s2;
double p1imm;
double p1phil;
double p2imm;
double p2phil;
double f1aimm;
double f1aphil;
double f1mimm;
double f1mphil;
double f2aimm;
double f2aphil;
double f2mimm;
double f2mphil;
double v1aimm;
double v1aphil;
double v1mimm;
double v1mphil;
double v2aimm;
double v2aphil;
double v2mimm;
double v2mphil;

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
    //
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double p1imm = ((struct rparams *) params)->p1imm;
double p1phil = ((struct rparams *) params)->p1phil;
double p2imm = ((struct rparams *) params)->p2imm;
double p2phil = ((struct rparams *) params)->p2phil;
double f1aimm = ((struct rparams *) params)->f1aimm;
double f1aphil = ((struct rparams *) params)->f1aphil;
double f1mimm = ((struct rparams *) params)->f1mimm;
double f1mphil = ((struct rparams *) params)->f1mphil;
double f2aimm = ((struct rparams *) params)->f2aimm;
double f2aphil = ((struct rparams *) params)->f2aphil;
double f2mimm = ((struct rparams *) params)->f2mimm;
double f2mphil = ((struct rparams *) params)->f2mphil;
double v1aimm = ((struct rparams *) params)->v1aimm;
double v1aphil = ((struct rparams *) params)->v1aphil;
double v1mimm = ((struct rparams *) params)->v1mimm;
double v1mphil = ((struct rparams *) params)->v1mphil;
double v2aimm = ((struct rparams *) params)->v2aimm;
double v2aphil = ((struct rparams *) params)->v2aphil;
double v2mimm = ((struct rparams *) params)->v2mimm;
double v2mphil = ((struct rparams *) params)->v2mphil;

    
    
    // 
double dfaimm1dttplus1, dfaimm2dttplus1, dfmimm1dttplus1, dfmimm2dttplus1, dfaphil1dttplus1, dfaphil2dttplus1, dfmphil1dttplus1, dfmphil2dttplus1;
for (int iter=0; iter<1e08; ++iter) {
dfaimm1dttplus1 = bound01(1 - f1aphil - f1mimm - f1mphil - f2aimm - f2aphil - f2mimm - f2mphil);

dfaimm2dttplus1 = bound01(f2aimm + 0.01*(-(f2aimm*mua*(1 - d*(f1aimm*(1 - p1imm) + f2mimm*(1 - p1imm) + f1aphil*(1 - p1phil) + f2mphil*(1 - p1phil) + f1mimm*(1 - p2imm) + f2aimm*(1 - p2imm) + f1mphil*(1 - p2phil) + f2aphil*(1 - p2phil)))) + mum*(d*f2mimm*(f1aimm*(1 - p1imm) + f2mimm*(1 - p1imm) + f1aphil*(1 - p1phil) + f2mphil*(1 - p1phil) + f1mimm*(1 - p2imm) + f2aimm*(1 - p2imm) + f1mphil*(1 - p2phil) + f2aphil*(1 - p2phil)) + d*f2mphil*(f1aimm*(1 - p1imm) + f2mimm*(1 - p1imm) + f1aphil*(1 - p1phil) + f2mphil*(1 - p1phil) + f1mimm*(1 - p2imm) + f2aimm*(1 - p2imm) + f1mphil*(1 - p2phil) + f2aphil*(1 - p2phil))) + d*f2aphil*mua*(f1aimm*(1 - p1imm) + f2mimm*(1 - p1imm) + f1aphil*(1 - p1phil) + f2mphil*(1 - p1phil) + f1mimm*(1 - p2imm) + f2aimm*(1 - p2imm) + f1mphil*(1 - p2phil) + f2aphil*(1 - p2phil)) + f1mimm*s1 - f2aimm*s2));

dfmimm1dttplus1 = bound01(f1mimm + 0.01*(-(f1mimm*mum*(1 - d*(f1aimm*(1 - p1imm) + f2mimm*(1 - p1imm) + f1aphil*(1 - p1phil) + f2mphil*(1 - p1phil) + f1mimm*(1 - p2imm) + f2aimm*(1 - p2imm) + f1mphil*(1 - p2phil) + f2aphil*(1 - p2phil)))) + mua*(d*f1aimm*(f1aimm*(1 - p1imm) + f2mimm*(1 - p1imm) + f1aphil*(1 - p1phil) + f2mphil*(1 - p1phil) + f1mimm*(1 - p2imm) + f2aimm*(1 - p2imm) + f1mphil*(1 - p2phil) + f2aphil*(1 - p2phil)) + d*f1aphil*(f1aimm*(1 - p1imm) + f2mimm*(1 - p1imm) + f1aphil*(1 - p1phil) + f2mphil*(1 - p1phil) + f1mimm*(1 - p2imm) + f2aimm*(1 - p2imm) + f1mphil*(1 - p2phil) + f2aphil*(1 - p2phil))) + d*f1mphil*mum*(f1aimm*(1 - p1imm) + f2mimm*(1 - p1imm) + f1aphil*(1 - p1phil) + f2mphil*(1 - p1phil) + f1mimm*(1 - p2imm) + f2aimm*(1 - p2imm) + f1mphil*(1 - p2phil) + f2aphil*(1 - p2phil)) - f1mimm*s1 + f2aimm*s2));

dfmimm2dttplus1 = bound01(f2mimm + 0.01*(d*f2mphil*mum*(f1aimm*p1imm + f2mimm*p1imm + f1aphil*p1phil + f2mphil*p1phil + f1mimm*p2imm + f2aimm*p2imm + f1mphil*p2phil + f2aphil*p2phil) - f2mimm*mum*(1 - d*(f1aimm*p1imm + f2mimm*p1imm + f1aphil*p1phil + f2mphil*p1phil + f1mimm*p2imm + f2aimm*p2imm + f1mphil*p2phil + f2aphil*p2phil)) + mua*(d*f2aimm*(f1aimm*p1imm + f2mimm*p1imm + f1aphil*p1phil + f2mphil*p1phil + f1mimm*p2imm + f2aimm*p2imm + f1mphil*p2phil + f2aphil*p2phil) + d*f2aphil*(f1aimm*p1imm + f2mimm*p1imm + f1aphil*p1phil + f2mphil*p1phil + f1mimm*p2imm + f2aimm*p2imm + f1mphil*p2phil + f2aphil*p2phil)) + f1aimm*s1 - f2mimm*s2));

dfaphil1dttplus1 = bound01(f1aphil + 0.01*((1 - d)*f1aimm*mua*p1imm - f1aphil*mua*(1 - (1 - d)*p1phil) + mum*((1 - d)*f1mimm*p2imm + (1 - d)*f1mphil*p2phil) - f1aphil*s1 + f2mphil*s2));

dfaphil2dttplus1 = bound01(f2aphil + 0.01*(mum*((1 - d)*f2mimm*(1 - p1imm) + (1 - d)*f2mphil*(1 - p1phil)) + (1 - d)*f2aimm*mua*(1 - p2imm) - f2aphil*mua*(1 - (1 - d)*(1 - p2phil)) + f1mphil*s1 - f2aphil*s2));

dfmphil1dttplus1 = bound01(f1mphil + 0.01*(mua*((1 - d)*f1aimm*(1 - p1imm) + (1 - d)*f1aphil*(1 - p1phil)) + (1 - d)*f1mimm*mum*(1 - p2imm) - f1mphil*mum*(1 - (1 - d)*(1 - p2phil)) - f1mphil*s1 + f2aphil*s2));

dfmphil2dttplus1 = bound01(f2mphil + 0.01*((1 - d)*f2mimm*mum*p1imm - f2mphil*mum*(1 - (1 - d)*p1phil) + mua*((1 - d)*f2aimm*p2imm + (1 - d)*f2aphil*p2phil) + f1aphil*s1 - f2mphil*s2));

if (
fabs(dfaimm1dttplus1 - f1aimm) < 1e-07 
&& fabs(dfaimm2dttplus1 - f2aimm) < 1e-07 
&& fabs(dfmimm1dttplus1 - f1mimm) < 1e-07 
&& fabs(dfmimm2dttplus1 - f2mimm) < 1e-07
&& fabs(dfaphil1dttplus1 - f1aphil) < 1e-07 
&& fabs(dfaphil2dttplus1 - f2aphil) < 1e-07 
&& fabs(dfmphil1dttplus1 - f1mphil) < 1e-07 
&& fabs(dfmphil2dttplus1 - f2mphil) < 1e-07) {
break;
}

f1aimm = dfaimm1dttplus1;
f2aimm = dfaimm2dttplus1;
f1mimm = dfmimm1dttplus1;
f2mimm = dfmimm2dttplus1;
f1aphil = dfaphil1dttplus1;
f2aphil = dfaphil2dttplus1;
f1mphil = dfmphil1dttplus1;
f2mphil = dfmphil2dttplus1;}

gsl_vector_set(f, 0, f1aimm);
gsl_vector_set(f, 1, f2aimm);
gsl_vector_set(f, 2, f1mimm);
gsl_vector_set(f, 3, f2mimm);
gsl_vector_set(f, 4, f1aphil);
gsl_vector_set(f, 5, f2aphil);
gsl_vector_set(f, 6, f1mphil);
gsl_vector_set(f, 7, f2mphil);
    
    return GSL_SUCCESS;
}


void reproductive_values(void *params, gsl_vector *f)
{
    // 
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double p1imm = ((struct rparams *) params)->p1imm;
double p1phil = ((struct rparams *) params)->p1phil;
double p2imm = ((struct rparams *) params)->p2imm;
double p2phil = ((struct rparams *) params)->p2phil;
double f1aimm = ((struct rparams *) params)->f1aimm;
double f1aphil = ((struct rparams *) params)->f1aphil;
double f1mimm = ((struct rparams *) params)->f1mimm;
double f1mphil = ((struct rparams *) params)->f1mphil;
double f2aimm = ((struct rparams *) params)->f2aimm;
double f2aphil = ((struct rparams *) params)->f2aphil;
double f2mimm = ((struct rparams *) params)->f2mimm;
double f2mphil = ((struct rparams *) params)->f2mphil;
double v1aimm = ((struct rparams *) params)->v1aimm;
double v1aphil = ((struct rparams *) params)->v1aphil;
double v1mimm = ((struct rparams *) params)->v1mimm;
double v1mphil = ((struct rparams *) params)->v1mphil;
double v2aimm = ((struct rparams *) params)->v2aimm;
double v2aphil = ((struct rparams *) params)->v2aphil;
double v2mimm = ((struct rparams *) params)->v2mimm;
double v2mphil = ((struct rparams *) params)->v2mphil;



    // 
double dvaimm1dttplus1, dvaimm2dttplus1, dvmimm1dttplus1, dvmimm2dttplus1, dvaphil1dttplus1, dvaphil2dttplus1, dvmphil1dttplus1, dvmphil2dttplus1;
for (int iter=0; iter<1e08; ++iter) {
dvaimm1dttplus1 = bound0(v1aimm + 0.01*(-(d*mua*v1aimm) + (1 - d)*mua*(-v1aimm + p1imm*v1aphil + (1 - p1imm)*v1mphil) + s1*(-v1aimm + v2mimm) + d*(f1aimm*mua*(p1imm*v1aimm + (1 - p1imm)*v1mimm) + f1aphil*mua*(p1imm*v1aimm + (1 - p1imm)*v1mimm) + f1mimm*mum*(p1imm*v1aimm + (1 - p1imm)*v1mimm) + f1mphil*mum*(p1imm*v1aimm + (1 - p1imm)*v1mimm) + f2aimm*mua*((1 - p1imm)*v2aimm + p1imm*v2mimm) + f2aphil*mua*((1 - p1imm)*v2aimm + p1imm*v2mimm) + f2mimm*mum*((1 - p1imm)*v2aimm + p1imm*v2mimm) + f2mphil*mum*((1 - p1imm)*v2aimm + p1imm*v2mimm))));

dvaimm2dttplus1 = bound0(v2aimm + 0.01*(s2*(v1mimm - v2aimm) - d*mua*v2aimm + d*(f1aimm*mua*(p2imm*v1aimm + (1 - p2imm)*v1mimm) + f1aphil*mua*(p2imm*v1aimm + (1 - p2imm)*v1mimm) + f1mimm*mum*(p2imm*v1aimm + (1 - p2imm)*v1mimm) + f1mphil*mum*(p2imm*v1aimm + (1 - p2imm)*v1mimm) + f2aimm*mua*((1 - p2imm)*v2aimm + p2imm*v2mimm) + f2aphil*mua*((1 - p2imm)*v2aimm + p2imm*v2mimm) + f2mimm*mum*((1 - p2imm)*v2aimm + p2imm*v2mimm) + f2mphil*mum*((1 - p2imm)*v2aimm + p2imm*v2mimm)) + (1 - d)*mua*(-v2aimm + (1 - p2imm)*v2aphil + p2imm*v2mphil)));

dvmimm1dttplus1 = bound0(v1mimm + 0.01*(-(d*mum*v1mimm) + (1 - d)*mum*(p2imm*v1aphil - v1mimm + (1 - p2imm)*v1mphil) + s1*(-v1mimm + v2aimm) + d*(f1aimm*mua*(p2imm*v1aimm + (1 - p2imm)*v1mimm) + f1aphil*mua*(p2imm*v1aimm + (1 - p2imm)*v1mimm) + f1mimm*mum*(p2imm*v1aimm + (1 - p2imm)*v1mimm) + f1mphil*mum*(p2imm*v1aimm + (1 - p2imm)*v1mimm) + f2aimm*mua*((1 - p2imm)*v2aimm + p2imm*v2mimm) + f2aphil*mua*((1 - p2imm)*v2aimm + p2imm*v2mimm) + f2mimm*mum*((1 - p2imm)*v2aimm + p2imm*v2mimm) + f2mphil*mum*((1 - p2imm)*v2aimm + p2imm*v2mimm))));

dvmimm2dttplus1 = bound0(v2mimm + 0.01*(s2*(v1aimm - v2mimm) - d*mum*v2mimm + d*(f1aimm*mua*(p1imm*v1aimm + (1 - p1imm)*v1mimm) + f1aphil*mua*(p1imm*v1aimm + (1 - p1imm)*v1mimm) + f1mimm*mum*(p1imm*v1aimm + (1 - p1imm)*v1mimm) + f1mphil*mum*(p1imm*v1aimm + (1 - p1imm)*v1mimm) + f2aimm*mua*((1 - p1imm)*v2aimm + p1imm*v2mimm) + f2aphil*mua*((1 - p1imm)*v2aimm + p1imm*v2mimm) + f2mimm*mum*((1 - p1imm)*v2aimm + p1imm*v2mimm) + f2mphil*mum*((1 - p1imm)*v2aimm + p1imm*v2mimm)) + (1 - d)*mum*((1 - p1imm)*v2aphil - v2mimm + p1imm*v2mphil)));

dvaphil1dttplus1 = bound0(v1aphil + 0.01*(-(d*mua*v1aphil) + (1 - d)*mua*(-v1aphil + p1phil*v1aphil + (1 - p1phil)*v1mphil) + d*(f1aimm*mua*(p1phil*v1aimm + (1 - p1phil)*v1mimm) + f1aphil*mua*(p1phil*v1aimm + (1 - p1phil)*v1mimm) + f1mimm*mum*(p1phil*v1aimm + (1 - p1phil)*v1mimm) + f1mphil*mum*(p1phil*v1aimm + (1 - p1phil)*v1mimm) + f2aimm*mua*((1 - p1phil)*v2aimm + p1phil*v2mimm) + f2aphil*mua*((1 - p1phil)*v2aimm + p1phil*v2mimm) + f2mimm*mum*((1 - p1phil)*v2aimm + p1phil*v2mimm) + f2mphil*mum*((1 - p1phil)*v2aimm + p1phil*v2mimm)) + s1*(-v1aphil + v2mphil)));

dvaphil2dttplus1 = bound0(v2aphil + 0.01*(s2*(v1mphil - v2aphil) - d*mua*v2aphil + d*(f1aimm*mua*(p2phil*v1aimm + (1 - p2phil)*v1mimm) + f1aphil*mua*(p2phil*v1aimm + (1 - p2phil)*v1mimm) + f1mimm*mum*(p2phil*v1aimm + (1 - p2phil)*v1mimm) + f1mphil*mum*(p2phil*v1aimm + (1 - p2phil)*v1mimm) + f2aimm*mua*((1 - p2phil)*v2aimm + p2phil*v2mimm) + f2aphil*mua*((1 - p2phil)*v2aimm + p2phil*v2mimm) + f2mimm*mum*((1 - p2phil)*v2aimm + p2phil*v2mimm) + f2mphil*mum*((1 - p2phil)*v2aimm + p2phil*v2mimm)) + (1 - d)*mua*(-v2aphil + (1 - p2phil)*v2aphil + p2phil*v2mphil)));

dvmphil1dttplus1 = bound0(v1mphil + 0.01*(-(d*mum*v1mphil) + (1 - d)*mum*(p2phil*v1aphil - v1mphil + (1 - p2phil)*v1mphil) + s1*(-v1mphil + v2aphil) + d*(f1aimm*mua*(p2phil*v1aimm + (1 - p2phil)*v1mimm) + f1aphil*mua*(p2phil*v1aimm + (1 - p2phil)*v1mimm) + f1mimm*mum*(p2phil*v1aimm + (1 - p2phil)*v1mimm) + f1mphil*mum*(p2phil*v1aimm + (1 - p2phil)*v1mimm) + f2aimm*mua*((1 - p2phil)*v2aimm + p2phil*v2mimm) + f2aphil*mua*((1 - p2phil)*v2aimm + p2phil*v2mimm) + f2mimm*mum*((1 - p2phil)*v2aimm + p2phil*v2mimm) + f2mphil*mum*((1 - p2phil)*v2aimm + p2phil*v2mimm))));

dvmphil2dttplus1 = bound0(v2mphil + 0.01*(d*(f1aimm*mua*(p1phil*v1aimm + (1 - p1phil)*v1mimm) + f1aphil*mua*(p1phil*v1aimm + (1 - p1phil)*v1mimm) + f1mimm*mum*(p1phil*v1aimm + (1 - p1phil)*v1mimm) + f1mphil*mum*(p1phil*v1aimm + (1 - p1phil)*v1mimm) + f2aimm*mua*((1 - p1phil)*v2aimm + p1phil*v2mimm) + f2aphil*mua*((1 - p1phil)*v2aimm + p1phil*v2mimm) + f2mimm*mum*((1 - p1phil)*v2aimm + p1phil*v2mimm) + f2mphil*mum*((1 - p1phil)*v2aimm + p1phil*v2mimm)) + s2*(v1aphil - v2mphil) - d*mum*v2mphil + (1 - d)*mum*((1 - p1phil)*v2aphil - v2mphil + p1phil*v2mphil)));

if (
fabs(dvaimm1dttplus1 - v1aimm) < 1e-07 
&& fabs(dvaimm2dttplus1 - v2aimm) < 1e-07 
&& fabs(dvmimm1dttplus1 - v1mimm) < 1e-07 
&& fabs(dvmimm2dttplus1 - v2mimm) < 1e-07
&& fabs(dvaphil1dttplus1 - v1aphil) < 1e-07 
&& fabs(dvaphil2dttplus1 - v2aphil) < 1e-07 
&& fabs(dvmphil1dttplus1 - v1mphil) < 1e-07 
&& fabs(dvmphil2dttplus1 - v2mphil) < 1e-07) {
break;
}

v1aimm = dvaimm1dttplus1;
v2aimm = dvaimm2dttplus1;
v1mimm = dvmimm1dttplus1;
v2mimm = dvmimm2dttplus1;
v1aphil = dvaphil1dttplus1;
v2aphil = dvaphil2dttplus1;
v1mphil = dvmphil1dttplus1;
v2mphil = dvmphil2dttplus1;}

gsl_vector_set(f, 0, v1aimm);
gsl_vector_set(f, 1, v2aimm);
gsl_vector_set(f, 2, v1mimm);
gsl_vector_set(f, 3, v2mimm);
gsl_vector_set(f, 4, v1aphil);
gsl_vector_set(f, 5, v2aphil);
gsl_vector_set(f, 6, v1mphil);
gsl_vector_set(f, 7, v2mphil);
}

void selgrads(void *params, gsl_vector *f)
{
    // 
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double p1imm = ((struct rparams *) params)->p1imm;
double p1phil = ((struct rparams *) params)->p1phil;
double p2imm = ((struct rparams *) params)->p2imm;
double p2phil = ((struct rparams *) params)->p2phil;
double f1aimm = ((struct rparams *) params)->f1aimm;
double f1aphil = ((struct rparams *) params)->f1aphil;
double f1mimm = ((struct rparams *) params)->f1mimm;
double f1mphil = ((struct rparams *) params)->f1mphil;
double f2aimm = ((struct rparams *) params)->f2aimm;
double f2aphil = ((struct rparams *) params)->f2aphil;
double f2mimm = ((struct rparams *) params)->f2mimm;
double f2mphil = ((struct rparams *) params)->f2mphil;
double v1aimm = ((struct rparams *) params)->v1aimm;
double v1aphil = ((struct rparams *) params)->v1aphil;
double v1mimm = ((struct rparams *) params)->v1mimm;
double v1mphil = ((struct rparams *) params)->v1mphil;
double v2aimm = ((struct rparams *) params)->v2aimm;
double v2aphil = ((struct rparams *) params)->v2aphil;
double v2mimm = ((struct rparams *) params)->v2mimm;
double v2mphil = ((struct rparams *) params)->v2mphil;


    // 
double dp1imm = bound01(p1imm + 0.01*(f1aimm*((1 - d)*mua*(v1aphil - v1mphil) + d*(f1aimm*mua*(v1aimm - v1mimm) + f1aphil*mua*(v1aimm - v1mimm) + f1mimm*mum*(v1aimm - v1mimm) + f1mphil*mum*(v1aimm - v1mimm) + f2aimm*mua*(-v2aimm + v2mimm) + f2aphil*mua*(-v2aimm + v2mimm) + f2mimm*mum*(-v2aimm + v2mimm) + f2mphil*mum*(-v2aimm + v2mimm))) + f2mimm*(d*(f1aimm*mua*(v1aimm - v1mimm) + f1aphil*mua*(v1aimm - v1mimm) + f1mimm*mum*(v1aimm - v1mimm) + f1mphil*mum*(v1aimm - v1mimm) + f2aimm*mua*(-v2aimm + v2mimm) + f2aphil*mua*(-v2aimm + v2mimm) + f2mimm*mum*(-v2aimm + v2mimm) + f2mphil*mum*(-v2aimm + v2mimm)) + (1 - d)*mum*(-v2aphil + v2mphil))));

double dp1phil = bound01(p1phil + 0.01*(f1aphil*((1 - d)*mua*(v1aphil - v1mphil) + d*(f1aimm*mua*(v1aimm - v1mimm) + f1aphil*mua*(v1aimm - v1mimm) + f1mimm*mum*(v1aimm - v1mimm) + f1mphil*mum*(v1aimm - v1mimm) + f2aimm*mua*(-v2aimm + v2mimm) + f2aphil*mua*(-v2aimm + v2mimm) + f2mimm*mum*(-v2aimm + v2mimm) + f2mphil*mum*(-v2aimm + v2mimm))) + f2mphil*(d*(f1aimm*mua*(v1aimm - v1mimm) + f1aphil*mua*(v1aimm - v1mimm) + f1mimm*mum*(v1aimm - v1mimm) + f1mphil*mum*(v1aimm - v1mimm) + f2aimm*mua*(-v2aimm + v2mimm) + f2aphil*mua*(-v2aimm + v2mimm) + f2mimm*mum*(-v2aimm + v2mimm) + f2mphil*mum*(-v2aimm + v2mimm)) + (1 - d)*mum*(-v2aphil + v2mphil))));

double dp2imm = bound01(p2imm + 0.01*(f1mimm*((1 - d)*mum*(v1aphil - v1mphil) + d*(f1aimm*mua*(v1aimm - v1mimm) + f1aphil*mua*(v1aimm - v1mimm) + f1mimm*mum*(v1aimm - v1mimm) + f1mphil*mum*(v1aimm - v1mimm) + f2aimm*mua*(-v2aimm + v2mimm) + f2aphil*mua*(-v2aimm + v2mimm) + f2mimm*mum*(-v2aimm + v2mimm) + f2mphil*mum*(-v2aimm + v2mimm))) + f2aimm*(d*(f1aimm*mua*(v1aimm - v1mimm) + f1aphil*mua*(v1aimm - v1mimm) + f1mimm*mum*(v1aimm - v1mimm) + f1mphil*mum*(v1aimm - v1mimm) + f2aimm*mua*(-v2aimm + v2mimm) + f2aphil*mua*(-v2aimm + v2mimm) + f2mimm*mum*(-v2aimm + v2mimm) + f2mphil*mum*(-v2aimm + v2mimm)) + (1 - d)*mua*(-v2aphil + v2mphil))));

double dp2phil = bound01(p2phil + 0.01*(f1mphil*((1 - d)*mum*(v1aphil - v1mphil) + d*(f1aimm*mua*(v1aimm - v1mimm) + f1aphil*mua*(v1aimm - v1mimm) + f1mimm*mum*(v1aimm - v1mimm) + f1mphil*mum*(v1aimm - v1mimm) + f2aimm*mua*(-v2aimm + v2mimm) + f2aphil*mua*(-v2aimm + v2mimm) + f2mimm*mum*(-v2aimm + v2mimm) + f2mphil*mum*(-v2aimm + v2mimm))) + f2aphil*(d*(f1aimm*mua*(v1aimm - v1mimm) + f1aphil*mua*(v1aimm - v1mimm) + f1mimm*mum*(v1aimm - v1mimm) + f1mphil*mum*(v1aimm - v1mimm) + f2aimm*mua*(-v2aimm + v2mimm) + f2aphil*mua*(-v2aimm + v2mimm) + f2mimm*mum*(-v2aimm + v2mimm) + f2mphil*mum*(-v2aimm + v2mimm)) + (1 - d)*mua*(-v2aphil + v2mphil))));

gsl_vector_set(f, 0, dp1imm);
gsl_vector_set(f, 1, dp1phil);
gsl_vector_set(f, 2, dp2imm);
gsl_vector_set(f, 3, dp2phil);


}


void write_params(void *params)
{
    // 
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double p1imm = ((struct rparams *) params)->p1imm;
double p1phil = ((struct rparams *) params)->p1phil;
double p2imm = ((struct rparams *) params)->p2imm;
double p2phil = ((struct rparams *) params)->p2phil;
double f1aimm = ((struct rparams *) params)->f1aimm;
double f1aphil = ((struct rparams *) params)->f1aphil;
double f1mimm = ((struct rparams *) params)->f1mimm;
double f1mphil = ((struct rparams *) params)->f1mphil;
double f2aimm = ((struct rparams *) params)->f2aimm;
double f2aphil = ((struct rparams *) params)->f2aphil;
double f2mimm = ((struct rparams *) params)->f2mimm;
double f2mphil = ((struct rparams *) params)->f2mphil;
double v1aimm = ((struct rparams *) params)->v1aimm;
double v1aphil = ((struct rparams *) params)->v1aphil;
double v1mimm = ((struct rparams *) params)->v1mimm;
double v1mphil = ((struct rparams *) params)->v1mphil;
double v2aimm = ((struct rparams *) params)->v2aimm;
double v2aphil = ((struct rparams *) params)->v2aphil;
double v2mimm = ((struct rparams *) params)->v2mimm;
double v2mphil = ((struct rparams *) params)->v2mphil;



    // 
DataFile << endl << endl  << "mum;" << mum << endl
 << "mua;" << mua << endl
 << "d;" << d << endl
 << "s1;" << s1 << endl
 << "s2;" << s2 << endl
 << endl;
}


void write_data(void *params, int time)
{
    // 
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double p1imm = ((struct rparams *) params)->p1imm;
double p1phil = ((struct rparams *) params)->p1phil;
double p2imm = ((struct rparams *) params)->p2imm;
double p2phil = ((struct rparams *) params)->p2phil;
double f1aimm = ((struct rparams *) params)->f1aimm;
double f1aphil = ((struct rparams *) params)->f1aphil;
double f1mimm = ((struct rparams *) params)->f1mimm;
double f1mphil = ((struct rparams *) params)->f1mphil;
double f2aimm = ((struct rparams *) params)->f2aimm;
double f2aphil = ((struct rparams *) params)->f2aphil;
double f2mimm = ((struct rparams *) params)->f2mimm;
double f2mphil = ((struct rparams *) params)->f2mphil;
double v1aimm = ((struct rparams *) params)->v1aimm;
double v1aphil = ((struct rparams *) params)->v1aphil;
double v1mimm = ((struct rparams *) params)->v1mimm;
double v1mphil = ((struct rparams *) params)->v1mphil;
double v2aimm = ((struct rparams *) params)->v2aimm;
double v2aphil = ((struct rparams *) params)->v2aphil;
double v2mimm = ((struct rparams *) params)->v2mimm;
double v2mphil = ((struct rparams *) params)->v2mphil;



    if (time < 0)
    {
        // 
DataFile << "time;f2aimm;v1mimm;f2mimm;v1aimm;f1mphil;p1imm;f1aphil;p2phil;v2aphil;f2mphil;p2imm;v2aimm;f1mimm;v2mphil;v2mimm;f1aimm;v1mphil;p1phil;v1aphil;f2aphil;" << endl;
    }

    // 
DataFile << time << ";" << f2aimm << ";" << 
v1mimm << ";" << 
f2mimm << ";" << 
v1aimm << ";" << 
f1mphil << ";" << 
p1imm << ";" << 
f1aphil << ";" << 
p2phil << ";" << 
v2aphil << ";" << 
f2mphil << ";" << 
p2imm << ";" << 
v2aimm << ";" << 
f1mimm << ";" << 
v2mphil << ";" << 
v2mimm << ";" << 
f1aimm << ";" << 
v1mphil << ";" << 
p1phil << ";" << 
v1aphil << ";" << 
f2aphil << ";" << 
 endl;
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
    // 
		paramstruct.mum = atof(argv[2]);
		paramstruct.mua = atof(argv[3]);
		paramstruct.d = atof(argv[4]);
		paramstruct.s1 = atof(argv[5]);
		paramstruct.s2 = atof(argv[6]);
		paramstruct.p1imm = atof(argv[7]);
		paramstruct.p1phil = atof(argv[8]);
		paramstruct.p2imm = atof(argv[9]);
		paramstruct.p2phil = atof(argv[10]);
		paramstruct.f1aimm = atof(argv[11]);
		paramstruct.f1aphil = atof(argv[12]);
		paramstruct.f1mimm = atof(argv[13]);
		paramstruct.f1mphil = atof(argv[14]);
		paramstruct.f2aimm = atof(argv[15]);
		paramstruct.f2aphil = atof(argv[16]);
		paramstruct.f2mimm = atof(argv[17]);
		paramstruct.f2mphil = atof(argv[18]);
		paramstruct.v1aimm = atof(argv[19]);
		paramstruct.v1aphil = atof(argv[20]);
		paramstruct.v1mimm = atof(argv[21]);
		paramstruct.v1mphil = atof(argv[22]);
		paramstruct.v2aimm = atof(argv[23]);
		paramstruct.v2aphil = atof(argv[24]);
		paramstruct.v2mimm = atof(argv[25]);
		paramstruct.v2mphil = atof(argv[26]);

   

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

