#ifndef TWOPATCH_H
#define TWOPATCH_H

enum 

class TwoPatch
{
    private:
        double m[2];
        double p;
        double mu[2];
        double c[2];

        void calculate_uv();

    public:
        TwoPatch(int argc, char **argv);

}; // end class TwoPatch

#endif
