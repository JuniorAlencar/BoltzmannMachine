#ifndef WRITE_JSON_H
#define WRITE_JSON_H


#include <nlohmann/json.hpp>
#include "network.h"

void write_json(const std::string &filename,
                const Rede &bm,
                const std::vector<double> &av_s,
                const std::vector<double> &av_ss,
                const std::vector<double> &av_sss,
                const VecDoub &bm_av_s,
                const VecDoub &bm_av_ss,
                const VecDoub_IO &bm_av_sss);

#endif // WRITE_JSON_H