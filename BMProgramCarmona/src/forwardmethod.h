#ifndef FORWARDMETHOD_H
#define FORWARDMETHOD_H

#include "network.h"
#include "energy.h"

void exact_solution (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const double beta);
void exact_solution_comp (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, Doub &E, Doub &E2, const double beta);
void exact_solution_cap (Rede &r, Doub &E, Doub &E2, const double beta);
void exact_solution_triplet (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, VecDoub_IO &av_sss, const double beta);
void metropolis_cap (Rede &r, const int t_eq, const int t_step,  const int relx, const int rept, const double T, 
						Doub &E, Doub &E2);
void metropolis_cap_mag (Rede &r, const int t_eq, const int t_step,  const int relx, const int rept, const double T, 
						Doub &E, Doub &E2, double &M, double &M2);
void metropolis_triplet (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, VecDoub_IO &av_sss, const int &t_eq, const int &t_meas,  
				const int &t_rep, const int &n_rep, const double &beta);

#endif
