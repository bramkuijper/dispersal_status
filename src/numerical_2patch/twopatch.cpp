#include "twopatch.hpp"

TwoPatch::TwoPatch(int argc, char **argv) :
    m{0.0,0.0}
    ,p{0.0}
    ,mu{0.0,0.0}
    ,c{0.0,0.0}
{
}

void TwoPatch::calculate_uv()
{
    // specify resident transition matrix
    double A[4][4] = {
        // row 1 of the matrix: 
        // transitions towards
        // (z1,e1)
        {(1.0 - m[0] + p * m[0])*(1.0 - mu[0])
            ,p * m[1] * (1.0 - mu[0])
                ,(1.0 - m[0] + m[0] * p) * mu[1]
                ,m[1] * p * mu[1]}
        // row 2 of the matrix:
        // transitions towards
        // (z1,e2)
        ,{(m[0]


    gsl_matrix_view m = gsl_matrix_view_array(
}

