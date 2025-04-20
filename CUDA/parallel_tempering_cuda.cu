#include <curand_kernel.h>
#include "../scripts/include/parallel_tempering_kernels.cuh" // inclui as declarações

__global__ void setup_curand(curandState *states, unsigned long seed) {
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(seed, id, 0, &states[id]);
}

__global__ void metropolis_step(
    int* spins_all, double* h_all, double* J_all, int2* bonds_all,
    int n_spins, int n_bonds, double beta,
    curandState *rand_states, int n_replicas
) {
    int rep_id = blockIdx.x;
    if (rep_id >= n_replicas) return;

    int idx = threadIdx.x;
    if (idx >= n_spins) return;

    int spin_offset = rep_id * n_spins;
    int bond_offset = rep_id * n_bonds;

    int s_i = spins_all[spin_offset + idx];
    double dE = 2.0 * s_i * h_all[spin_offset + idx];

    for (int b = 0; b < n_bonds; ++b) {
        int2 bond = bonds_all[bond_offset + b];
        if (bond.x == idx || bond.y == idx) {
            int neighbor = (bond.x == idx) ? bond.y : bond.x;
            dE += 2.0 * J_all[bond_offset + b] * s_i * spins_all[spin_offset + neighbor];
        }
    }

    double r = curand_uniform(&rand_states[rep_id * n_spins + idx]);
    if (dE <= 0 || r < exp(-beta * dE)) {
        spins_all[spin_offset + idx] = -s_i;
    }
}

__global__ void accumulate_observables(
    int* spins_all, double* av_s_all, double* av_ss_all, int2* bonds_all,
    int n_spins, int n_bonds, int n_replicas
) {
    int rep_id = blockIdx.x;
    if (rep_id >= n_replicas) return;

    int idx = threadIdx.x;
    int spin_offset = rep_id * n_spins;
    int bond_offset = rep_id * n_bonds;

    if (idx < n_spins)
        atomicAdd(&av_s_all[spin_offset + idx], (double)spins_all[spin_offset + idx]);

    if (idx < n_bonds) {
        int2 bond = bonds_all[bond_offset + idx];
        int si = spins_all[spin_offset + bond.x];
        int sj = spins_all[spin_offset + bond.y];
        atomicAdd(&av_ss_all[bond_offset + idx], (double)(si * sj));
    }
}

__global__ void compute_energies_kernel(
    const int* spins_all, const double* h_all, const double* J_all, const int2* bonds_all,
    double* energies, int n_spins, int n_bonds, int n_replicas
) {
    int idx = blockIdx.x;
    if (idx >= n_replicas) return;

    int spin_offset = idx * n_spins;
    int bond_offset = idx * n_bonds;

    double E = 0.0;

    for (int i = 0; i < n_spins; ++i)
        E -= h_all[spin_offset + i] * spins_all[spin_offset + i];

    for (int i = 0; i < n_bonds; ++i) {
        int2 b = bonds_all[bond_offset + i];
        int s1 = spins_all[spin_offset + b.x];
        int s2 = spins_all[spin_offset + b.y];
        E -= J_all[bond_offset + i] * s1 * s2;
    }

    energies[idx] = E;
}




