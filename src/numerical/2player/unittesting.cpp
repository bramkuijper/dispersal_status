#define CATCH_CONFIG_MAIN 
#include <iostream>
#include <catch2/catch.hpp>
#include "patchfreq.hpp"

// tests for the components of our programme


TEST_CASE("Check basic functionality of patch frequency ordered arrays","[patchfreq]") {

    PatchFreq pf{};

    for (int envt_i = 0; envt_i < 2; ++envt_i)
    {
        for (int state1_i = 0; state1_i < 4; ++state1_i)
        {
            for (int state2_i = 0; state2_i < 4; ++state2_i)
            {
                INFO("The state is " << envt_i << " " << state1_i << " " << state2_i);
                CHECK(pf(envt_i, state1_i, state2_i) == 1.0/20);
            }
        }
    }

    // set one of the frequencies to a higher value
    pf(0,1,2) = 0.5;
    pf(1,2,3) = 0.8;

    CHECK(pf(0,1,2) == pf(0,2,1));
    CHECK(pf(1,2,3) == pf(1,2,3));
}
