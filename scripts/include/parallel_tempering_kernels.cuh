#ifndef PARALLEL_TEMPERING_KERNELS_CUH
#define PARALLEL_TEMPERING_KERNELS_CUH

#include <curand_kernel.h>


__global__ void metropolis_step(int* spins_all, double* h_all, double* J_all, int2* bonds_all,
                                int n_spins, int n_bonds, double beta,
                                curandState *rand_states, int n_replicas);

__global__ void accumulate_observables(int* spins_all, double* av_s_all, double* av_ss_all, int2* bonds_all,
                                       int n_spins, int n_bonds, int n_replicas);

__global__ void compute_energies_kernel(const int* spins_all, const double* h_all, const double* J_all, const int2* bonds_all,
                                        double* energies, int n_spins, int n_bonds, int n_replicas);

__global__ void setup_curand(curandState *states, unsigned long seed);

#endif
