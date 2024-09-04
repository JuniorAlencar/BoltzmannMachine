#ifndef WRITE_JSON_H
#define WRITE_JSON_H


#include <nlohmann/json.hpp>
#include "network.h"

void write_json(const std::string &filename,
                const Rede &bm, // Network
                const std::vector<double> &av_s, // first moment exp
                const std::vector<double> &av_ss, // second moment exp
                const std::vector<double> &av_sss, // third moment exp
                const VecDoub &bm_av_s, // first moment ising
                const VecDoub &bm_av_ss, // second moment ising
                const VecDoub_IO &bm_av_sss); // third moment ising

// void write_properties(const std::string &filename,
//                        const Rede &bm
// );
// void write_specific_heat(const std::string &filename,)
//                         ;

#endif // WRITE_JSON_H