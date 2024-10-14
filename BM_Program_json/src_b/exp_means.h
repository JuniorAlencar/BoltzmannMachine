#ifndef EXP_MEANS_H
#define EXP_MEANS_H

#include "json_functions.h"

class exp_mean_calculate{
    private:
        bool use_exact;             // method exact/metropolis
        string filename;
    
    public:
        exp_means exp_calculate(const string &filename);
};


#endif // !EXP_MEANS_H
