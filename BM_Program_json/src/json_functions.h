#ifndef WRITE_JSON_H
#define WRITE_JSON_H

#include "create_folders.h"
#include <iostream>
#include <nlohmann/json.hpp> // Inclua a biblioteca JSON antes das outras
#include "nr3.h"
#include "network.h"

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

class js_funct{
    public:
        void create_json_exp(const exp_means &data, const string &file_means);
        void write_json_properties(const string &filename,                 // Sample name
                            const int &nspins,
                            const int &n_duplet,                         // Network
                            // MC parameters
                            const int &multi_relx,                        // Parameter to MC
                            const int &multi_teq,                         // Parameter to MC
                            const double &jmin,                     // Parameter to MC
                            const double &hmin,                     // Parameter to MC
                            const string &method,                   // exact or metropolis
                            // Exp properties
                            const exp_means &data_exp,              //experimental means
                            // Ising properties
                            const VecDoub &bm_av_s,                 // First moment Ising
                            const VecDoub &bm_av_ss,                // Second moment Ising
                            const std::vector<double> &bm_av_sss,   // Third moment Ising
                            const std::vector<double> &Pij_ising,   // Correlation Ising
                            const std::vector<double> &Tijk_ising   // Triplet Ising
                            );
        exp_means load_json_exp(const std::string &file_means);
};

#endif // !WRITE_JSON_H
