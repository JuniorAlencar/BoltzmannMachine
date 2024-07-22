#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "network.h"
#include "energy.h"

void metropolis (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, VecDoub_IO &av_sss, const int &t_eq, const int &t_meas,  
					const int &t_rep, const int &n_rep, const double &beta, const bool &tr);

#endif //METROPOLIS_H