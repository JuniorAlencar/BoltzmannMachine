#include <cuda_runtime.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>
#include <numeric>
#include <cuda_runtime.h>
#include "../scripts/include/cuda_tools.h"
#include "../scripts/include/nr3.h"
#include "../scripts/include/network.h"
#include "../scripts/include/forwardmethod_jr.h"
#include "../scripts/include/parallel_tempering_kernels.cuh"

void parallel_tempering_cuda_multi(
    int n_replicas, double T_min, double T_max,
    VecDoub_IO &bm_av_s, VecDoub_IO &bm_av_ss,
    int t_eq, int t_step, int relx, int rept,
    int n_spins, double mean, double sigma,
    const int &type, double H,
    std::vector<double> &energy_per_replica,
    std::vector<double> &temperatures,
    double &swap_acceptance_ratio,
    std::mt19937 &gen
) {
    temperatures.resize(n_replicas);
    int idx_T1 = n_replicas / 2;
    temperatures[idx_T1] = 1.0;
    for (int i = 0; i < idx_T1; ++i)
        temperatures[i] = 1.0 * pow(T_min / 1.0, (double)(idx_T1 - i) / idx_T1);
    for (int i = idx_T1 + 1; i < n_replicas; ++i)
        temperatures[i] = 1.0 * pow(T_max / 1.0, (double)(i - idx_T1) / (n_replicas - 1 - idx_T1));

    std::vector<double> betas(n_replicas);
    std::vector<Rede> replicas;
    replicas.reserve(n_replicas);
    for (int i = 0; i < n_replicas; ++i) {
        betas[i] = 1.0 / temperatures[i];
        replicas.emplace_back(n_spins, mean, sigma, betas[i], type, H);
        replicas.back().create_bonds_random();
    }

    int n_bonds = replicas[0].get_all_bonds().size() / 2;
    size_t sz_spins = n_replicas * n_spins * sizeof(int);
    size_t sz_h = n_replicas * n_spins * sizeof(double);
    size_t sz_J = n_replicas * n_bonds * sizeof(double);
    size_t sz_bonds = n_replicas * n_bonds * sizeof(int2);

    int *d_spins_all;
    double *d_h_all, *d_J_all, *d_energies;
    int2 *d_bonds_all;
    double *d_av_s, *d_av_ss;
    curandState *d_rand_states;

    cudaMalloc(&d_spins_all, sz_spins);
    cudaMalloc(&d_h_all, sz_h);
    cudaMalloc(&d_J_all, sz_J);
    cudaMalloc(&d_bonds_all, sz_bonds);
    cudaMalloc(&d_energies, sizeof(double) * n_replicas);
    cudaMalloc(&d_av_s, sz_spins);
    cudaMalloc(&d_av_ss, sz_J);
    cudaMalloc(&d_rand_states, n_replicas * n_spins * sizeof(curandState));

    std::vector<int> spins_all(n_replicas * n_spins);
    std::vector<double> h_all(n_replicas * n_spins);
    std::vector<double> J_all(n_replicas * n_bonds);
    std::vector<int2> bonds_all(n_replicas * n_bonds);

    for (int r = 0; r < n_replicas; ++r) {
        const std::vector<int>& bond_indices = replicas[r].get_all_bonds();
        for (int i = 0; i < n_spins; ++i) {
            spins_all[r * n_spins + i] = replicas[r].s[i];
            h_all[r * n_spins + i] = replicas[r].h[i];
        }
        for (int i = 0; i < n_bonds; ++i) {
            J_all[r * n_bonds + i] = replicas[r].J[i];
            int site1 = bond_indices[2 * i];
            int site2 = bond_indices[2 * i + 1];
            bonds_all[r * n_bonds + i] = make_int2(site1, site2);
        }
    }

    cudaMemcpy(d_spins_all, spins_all.data(), sz_spins, cudaMemcpyHostToDevice);
    cudaMemcpy(d_h_all, h_all.data(), sz_h, cudaMemcpyHostToDevice);
    cudaMemcpy(d_J_all, J_all.data(), sz_J, cudaMemcpyHostToDevice);
    cudaMemcpy(d_bonds_all, bonds_all.data(), sz_bonds, cudaMemcpyHostToDevice);
    cudaMemset(d_av_s, 0, sz_spins);
    cudaMemset(d_av_ss, 0, sz_J);

    setup_curand<<<n_replicas, n_spins>>>(d_rand_states, gen());

    std::vector<int> indices(n_replicas);
    std::iota(indices.begin(), indices.end(), 0);

    int swap_attempts = 0, swap_accepted = 0;

    for (int rep = 0; rep < rept; ++rep) {
        metropolis_step<<<n_replicas, n_spins>>>(
            d_spins_all, d_h_all, d_J_all, d_bonds_all,
            n_spins, n_bonds, 1.0, d_rand_states, n_replicas
        );
        cudaDeviceSynchronize();

        compute_energies_kernel<<<n_replicas, 1>>>(
            d_spins_all, d_h_all, d_J_all, d_bonds_all,
            d_energies, n_spins, n_bonds, n_replicas
        );
        cudaDeviceSynchronize();

        std::vector<double> energies(n_replicas);
        cudaMemcpy(energies.data(), d_energies, sizeof(double) * n_replicas, cudaMemcpyDeviceToHost);

        for (int i = 0; i < n_replicas - 1; ++i) {
            int idx_i = indices[i];
            int idx_j = indices[i + 1];
            double delta = (betas[idx_j] - betas[idx_i]) * (energies[idx_j] - energies[idx_i]);
            ++swap_attempts;
            if ((double)gen() / gen.max() < exp(delta)) {
                std::swap(indices[i], indices[i + 1]);
                ++swap_accepted;
            }
        }

        accumulate_observables<<<n_replicas, max(n_spins, n_bonds)>>>(
            d_spins_all, d_av_s, d_av_ss, d_bonds_all,
            n_spins, n_bonds, n_replicas
        );
        cudaDeviceSynchronize();
    }

    std::vector<double> av_s_host(n_replicas * n_spins);
    std::vector<double> av_ss_host(n_replicas * n_bonds);

    cudaMemcpy(av_s_host.data(), d_av_s, sz_spins, cudaMemcpyDeviceToHost);
    cudaMemcpy(av_ss_host.data(), d_av_ss, sz_J, cudaMemcpyDeviceToHost);

    int target_index = idx_T1;
    for (int i = 0; i < n_spins; ++i)
        bm_av_s[i] = av_s_host[target_index * n_spins + i] / (rept * t_step / relx);

    for (int i = 0; i < n_bonds; ++i)
        bm_av_ss[i] = av_ss_host[target_index * n_bonds + i] / (rept * t_step / relx);

    swap_acceptance_ratio = (swap_attempts > 0) ? (double)swap_accepted / swap_attempts : 0.0;

    cudaFree(d_spins_all);
    cudaFree(d_h_all);
    cudaFree(d_J_all);
    cudaFree(d_bonds_all);
    cudaFree(d_energies);
    cudaFree(d_av_s);
    cudaFree(d_av_ss);
    cudaFree(d_rand_states);
}





