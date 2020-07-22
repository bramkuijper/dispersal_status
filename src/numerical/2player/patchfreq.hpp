#ifndef PATCHFREQ_
#define PATCHFREQ_

enum State
{
    d1 = 0, // dispersal z1
    d2 = 1, // dispersal z2
    n1 = 2, // philopatric z1
    n2 = 3 // philopatric z2
};



class PatchFreq
{
    private:
        // size of the data structure
        static const int fsize = 20;

        // the data, an array of 10 values 
        // for the two player game
        double f[fsize];

    public:

        // initializes the arra
        PatchFreq();

        double & operator() (bool const envt_i, int const state1, int const state2);


};


#endif
