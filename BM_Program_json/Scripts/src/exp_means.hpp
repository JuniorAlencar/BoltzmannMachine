#ifndef EXP_MEANS_HPP
#define EXP_MEANS_HPP

#include <vector>
#include "create_folders.cpp"
#include "read_json.hpp"
#include "write_json.hpp"

using json = nlohmann::json;
using namespace std;

struct exp_means {
    std::vector<double> av_s;
    std::vector<double> av_ss;
    std::vector<double> av_sss;
    std::vector<double> Cij_exp;
    std::vector<double> Pij_exp;
    std::vector<double> Tijk_exp;
};

class exp_mean_calculate{
    private:
        bool use_exact;             // method exact/metropolis
        string filename;
    
    public:
        pair<exp_means, json> exp_calculate(string &filename, bool &use_exact);
};


#endif // !EXP_MEANS_HPP