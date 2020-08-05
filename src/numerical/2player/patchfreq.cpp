#include <iostream>
#include "patchfreq.hpp"

PatchFreq::PatchFreq()
{
    for (int f_idx = 0; f_idx < fsize; ++f_idx)
    {
        f[f_idx] = 1.0 / fsize;
    }
} 

// data get/set function
double & PatchFreq::operator() (
        bool const envt_i
        ,int const state1
        ,int const state2)
{
    // use ordering, so if state2 is larger than state1
    // then switch order
    if (state2 > state1)
    {
        assert(envt_i * 10 + state2 * (state2 + 1)/2 + state1 >= 0);
        assert(envt_i * 10 + state2 * (state2 + 1)/2 + state1 < fsize);
        // see Haldane 1948 J Genet. doi:10.1007/BF02986828 
        // Ordered combinations of genotypes
        // total number is 0.5 * g * (g+1)

        return(f[envt_i * 10 + state2 * (state2 + 1)/2 + state1]);
    }

    assert(envt_i * 10 + state1 * (state1 + 1)/2 + state2 >= 0);
    assert(envt_i * 10 + state1 * (state1 + 1)/2 + state2 < fsize);

//    std::cout << envt_i << " " << state1 << " " << state2 << " " << envt_i * 10 + (int) 0.5 * state1 * (state1 + 1) << " " << envt_i * 10 + state1 * (state1 + 1)/2 + state2 << std::endl;
    // state1 >= state2
    return(f[envt_i * 10 + state1 * (state1 + 1)/2 + state2]);
} // end elem

