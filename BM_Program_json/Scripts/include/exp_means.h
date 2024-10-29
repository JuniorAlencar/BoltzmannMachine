#ifndef EXP_MEANS_H
#define EXP_MEANS_H

#include <vector>

struct exp_means {
    std::vector<double> av_s;
    std::vector<double> av_ss;
    std::vector<double> av_sss;
    std::vector<double> Cij_exp;
    std::vector<double> Pij_exp;
    std::vector<double> Tijk_exp;
};

#endif // EXP_MEANS_H