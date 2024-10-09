#ifndef EXP_MEANS_HPP
#define EXP_MEANS_HPP

#include "json_functions.hpp"


using json = nlohmann::json;
using namespace std;

class exp_mean_calculate{
    private:
        bool use_exact;             // method exact/metropolis
        string filename;
    
    public:
        exp_means exp_calculate(const string &filename);
};


#endif // !EXP_MEANS_HPP