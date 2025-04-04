#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <bits/stdc++.h>
#include <cmath>
#include <numeric>
#include <optional>
#include "network.h"
#include "nr3.h"

using namespace std;

// If you use vscode, install the Doxygen Documentation Generator extension

// Compute properties =========================================================/

/**
 * @brief Generates a symmetric matrix J_ij with Gaussian-distributed values.
 *
 * @param N \c int: Number of spins
 * @param mean \c double: Mean of the Gaussian distribution
 * @param stddev \c double: Standard deviation of the Gaussian distribution
 * @param seed \c std::optional<unsigned int> (optional): Seed for the random number generator. If not provided, a random device is used.
 * @return Matrix J of size NxN with symmetric Gaussian-distributed values
 */
vector<vector<double>> generateJMatrix(int N, double mean = 0.0, double stddev, optional<unsigned int> seed = nullopt) {
    mt19937 gen(seed ? mt19937(*seed) : mt19937(random_device{}()));
    normal_distribution<double> dist(mean, stddev);

    vector<vector<double>> J(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            J[i][j] = J[j][i] = dist(gen);
        }
    }
    return J;
}


/**
 * @brief Generates a vector of external field values h_i from a uniform distribution.
 *
 * @param N \c int: Number of spins
 * @param min \c double: Minimum value of the external field
 * @param max \c double: Maximum value of the external field
 * @param seed \c std::optional<unsigned int>: Seed for the random number generator. If not provided, a random device is used.
 * @return Vector of size N with uniformly distributed values in [min, max]
 */
vector<double> generateHVector(int N, double min = -1.0, double max = 1.0, optional<unsigned int> seed = nullopt) {
    mt19937 gen(seed ? mt19937(*seed) : mt19937(random_device{}()));
    uniform_real_distribution<double> dist(min, max);

    vector<double> vec(N);
    for (int i = 0; i < N; ++i) {
        vec[i] = dist(gen);
    }
    return vec;
}

/**
 * @brief Runs the Monte Carlo simulation using the Metropolis algorithm for a spin model.
 *
 * @param r \c Rede: Structure containing spin vector, external field, and coupling matrix.
 * @param av_s \c VecDoub_IO: Accumulator for the average of individual spins over sampling.
 * @param av_ss \c VecDoub_IO: Accumulator for the average of pairwise spin products.
 * @param t_eq \c int: Number of steps for the system to reach equilibrium (relaxation time).
 * @param t_step \c int: Number of steps after equilibrium used for sampling.
 * @param relx \c int: Interval between spin state samplings after the system equilibrates.
 * @param rept \c int: Number of independent repetitions of the Monte Carlo process.
 * @param beta \c double: Inverse temperature \c (1/kT).
 * @param energies \c std::vector<double>: Vector storing the total energy of the system at each step.
 * @param states \c std::vector<std::vector<int>>: Matrix of size MxN storing the spin states over time.
 *               Each row corresponds to a sampled time step, and each column to a spin.
 * @param seed \c std::optional<unsigned int>: \c Optional seed for the random number generator.
 *              If provided, ensures reproducible simulations.
 */
void GenerateStates(
    Rede &r,
    VecDoub_IO &av_s,
    VecDoub_IO &av_ss,
    const int t_eq,
    const int t_step,
    const int relx,
    const int rept,
    const double beta,
    std::vector<double> &energies,
    std::vector<std::vector<int>> &states,
    std::optional<unsigned int> seed = std::nullopt
) {
    int N = r.n;
    int total_steps = t_eq + t_step;

    std::mt19937 gen(seed ? std::mt19937(*seed) : std::mt19937(std::random_device{}()));
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::uniform_int_distribution<int> spin_dist(0, N - 1);

    int M = t_step / relx;
    states.assign(M, std::vector<int>(N, 0)); // M rows = sampled states, N columns = spins
    energies.clear();
    int m_index = 0;

    for (int j = 0; j < rept; j++) {
        m_index = 0;

        for (int i = 0; i < total_steps; i++) {
            int k = spin_dist(gen);

            // Energy change if spin k is flipped
            double dE = 2.0 * r.s[k] * r.h[k];
            for (int l = 0; l < N; l++) {
                dE += 2.0 * r.s[k] * r.J[k][l] * r.s[l];
            }

            if (dE <= 0.0 || dist(gen) < std::exp(-beta * dE)) {
                r.s[k] *= -1;
            }

            // Compute total energy and store
            double energy = 0.0;
            for (int m = 0; m < N; m++) {
                energy -= r.h[m] * r.s[m];
                for (int n = m + 1; n < N; n++) {
                    energy -= r.J[m][n] * r.s[m] * r.s[n];
                }
            }
            energies.push_back(energy);

            // After equilibration, sample state every 'relx' steps
            if (i >= t_eq && ((i - t_eq) % relx == 0)) {
                for (int jj = 0; jj < N; jj++) {
                    av_s[jj] += r.s[jj];
                }

                int ind_ss = 0;
                for (int jj = 0; jj < N - 1; jj++) {
                    for (int l = jj + 1; l < N; l++) {
                        av_ss[ind_ss] += r.s[jj] * r.s[l];
                        ind_ss++;
                    }
                }

                // Store current spin state
                for (int spin = 0; spin < N; ++spin) {
                    states[m_index][spin] = r.s[spin];
                }
                m_index++;
            }
        }
    }
}

/**
 * @brief Computes the first moment (magnetization) of the system.
 *
 * @param sigma_states \c vector<vector<double>>: Matrix containing all spin states (M samples × N spins)
 * @return \c vector<double>: Magnetization vector s_i
 */
vector<double> computeSi(const vector<vector<double>>& sigma_states) {
    int N_samples = sigma_states.size(); // Number of rows (M_states)
    int N_spins = (N_samples > 0) ? sigma_states[0].size() : 0; // Number of columns (N_spins)
    
    vector<double> s(N_spins, 0.0);
    
    for (int p = 0; p < N_spins; p++) {
        for (int w = 0; w < N_samples; w++) {
            s[p] += sigma_states[w][p]; // Summing over all samples
        }
        s[p] /= N_samples;
    }

    return s;
}

/**
 * @brief Computes the second moment (⟨s_i s_j⟩) of the system.
 *
 * @param sigma_states \c vector<vector<double>>: Matrix containing all spin states (M samples × N spins)
 * @return \c vector<double>: Vector of s_i * s_j averages for all unique pairs (i < j)
 */
vector<double> computeSiSj(const vector<vector<double>>& sigma_states) {
    int N_samples = sigma_states.size();
    int N_spins = (N_samples > 0) ? sigma_states[0].size() : 0;
    int N_duplet = (N_spins * (N_spins - 1)) / 2;
    
    vector<double> ss(N_duplet, 0.0); // Pairwise product averages

    int aux = 0;
    for (int p = 0; p < N_spins - 1; p++) {
        for (int pp = p + 1; pp < N_spins; pp++) {
            for (int w = 0; w < N_samples; w++) {
                ss[aux] += sigma_states[w][p] * sigma_states[w][pp];
            }
            ss[aux] /= N_samples;
            aux++;
        }
    }

    return ss;
}

/**
 * @brief Computes the connected correlations C_ij = ⟨s_i s_j⟩ - ⟨s_i⟩⟨s_j⟩.
 *
 * @param s \c vector<double>: Magnetization vector ⟨s_i⟩
 * @param ss \c vector<double>: Second moment vector ⟨s_i s_j⟩
 * @return \c vector<double>: Connected correlation values C_ij
 */
vector<double> computeC(const vector<double> s, const vector<double> ss) {
    int N_spins = s.size();
    int N_duplet = ss.size();
    
    vector<double> C(N_duplet, 0.0);
    int aux = 0;
    for (int p = 0; p < N_spins - 1; p++) {
        for (int pp = p + 1; pp < N_spins; pp++) {
            C[aux] = ss[aux] - s[p] * s[pp];
            aux++;
        }
    }
    return C;
}

/**
 * @brief Computes the Hamiltonian H for each spin configuration in the system.
 *
 * @param sigmaStates \c vector<vector<double>>: Matrix of spin states (M samples × N spins)
 * @param h \c vector<double>: External field vector h_i
 * @param J \c vector<vector<double>>: Symmetric interaction matrix J_ij
 * @return \c vector<double>: Vector of Hamiltonian values for each configuration
 */
vector<double> computeHamiltonian(const vector<vector<double>>& sigmaStates, const vector<double>& h, const vector<vector<double>>& J) {
    int M = sigmaStates.size();
    int N = (M > 0) ? sigmaStates[0].size() : 0;
    vector<double> H_values(M, 0.0);
    
    for (int m = 0; m < M; ++m) {
        double H = 0.0;

        // External field contribution
        for (int i = 0; i < N; ++i) {
            H += h[i] * sigmaStates[m][i];
        }

        // Pairwise interaction term (only upper triangle considered)
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                H -= J[i][j] * sigmaStates[m][i] * sigmaStates[m][j];
            }
        }

        H_values[m] = H;
    }
    return H_values;
}


/**
 * @brief Computes magnetization (si), second moment (sisj), and correlation (C) from spin states.
 *
 * @param sigmaStates \c vector<vector<double>>: Matrix containing all the spin states (M samples × N spins)
 * @param si \c vector<double>&: Output vector of mean spin values (first moment)
 * @param sisj \c vector<double>&: Output vector of mean products s_i * s_j (second moment)
 * @param C \c vector<double>&: Output vector of connected correlations C_ij
 */
void computeMagCorr(const vector<vector<double>>& sigmaStates, vector<double>& si, vector<double>& sisj, vector<double>& C) {
    int M = sigmaStates.size();
    int N = sigmaStates[0].size();

    si.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
        for (int m = 0; m < M; ++m)
            si[i] += sigmaStates[m][i];
        si[i] /= M;
    }

    for (int i = 0; i < N - 1; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double mean_prod = 0.0;
            for (int m = 0; m < M; ++m)
                mean_prod += sigmaStates[m][i] * sigmaStates[m][j];
            mean_prod /= M;

            sisj.push_back(mean_prod);
            C.push_back(mean_prod - si[i] * si[j]);
        }
    }
}

// Save properties =========================================================/

/**
 * @brief Saves h_i values (external fields) to a file in column format.
 *
 * @param filename \c std::string: Path + filename (including extension)
 * @param data \c vector<double>: Vector containing h_i values
 */
void savehH(const string& filename, const vector<double>& data) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error opening file " << filename << endl;
        return;
    }

    for (const auto& value : data) {
        file << fixed << setprecision(4) << value << "\n";
    }
    file.close();
    cout << "File saved: " << filename << endl;
}

/**
 * @brief Saves the upper triangular part of the interaction matrix J_ij to file.
 *
 * @param filename \c  std::string: Path + filename (including extension)
 * @param J \c std::vector<vector<double>>: Symmetric interaction matrix J_ij
 */
void saveJ(const string& filename, const vector<vector<double>>& J) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error opening file " << filename << endl;
        return;
    }

    for (int i = 0; i < J.size(); ++i) {
        for (int j = i + 1; j < J[i].size(); ++j) { // Upper triangle only
            file << fixed << setprecision(4) << J[i][j];
            if (j < J[i].size() - 1) {
                file << ","; // Comma separator
            }
        }
        file << "\n"; // New row for each i-th row
    }

    file.close();
    cout << "File saved: " << filename << endl;
}

/**
 * @brief Saves multiple spin configurations (sigma states) row-wise in CSV-like format.
 *
 * @param filename \c std::string: Path + filename (including extension)
 * @param sigmaStates \c vector<vector<double>>: Matrix with spin states (each row = one configuration)
 */
void saveSigmaStates(const string& filename, const vector<vector<double>>& sigmaStates) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error opening file " << filename << endl;
        return;
    }

    for (const auto& state : sigmaStates) {
        for (size_t i = 0; i < state.size(); ++i) {
            file << fixed << setprecision(1) << state[i]; // One decimal place
            if (i < state.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
    cout << "File saved: " << filename << endl;
}

/**
 * @brief Saves the magnetization vector si (first moment) to a file.
 *
 * @param filename \c std::string: Name of file, including path and extension
 * @param si \c std::vector<double>: Vector containing the mean value of each spin (s_i)
 */
void saveSi(const string& filename, const vector<double>& si) {
    ofstream file_si(filename.c_str());
    int N_spins = si.size();

    for (int i = 0; i < N_spins; i++) {
        file_si << si[i] << endl;
    }

    file_si.close();
    cout << "File saved: " << filename << endl;
}

/**
 * @brief Saves the second moment vector sisj (⟨s_i s_j⟩) to a file.
 *
 * @param filename \c std::string: Name of file, including path and extension
 * @param sisj \c vector<double>: Vector of pairwise spin products
 */
void saveSiSj(const string& filename, const vector<double>& sisj) {
    ofstream file_sisj(filename.c_str());
    int N_duplet = sisj.size();

    for (int i = 0; i < N_duplet; i++) {
        file_sisj << sisj[i] << endl;
    }

    file_sisj.close();
    cout << "File saved: " << filename << endl;
}

/**
 * @brief Saves si, sisj, and C vectors to a file in a compatible format.
 *
 * The file starts with the number of spins, followed by one line per duplet:
 *   <sisj> <C> <si> (if si index is within bounds)
 *
 * @param filename \c string: Name of file, including path and extension
 * @param si \c vector<double>: First moment (mean values of spins)
 * @param sisj \c vector<double>: Second moment (⟨s_i s_j⟩)
 * @param C \c vector<double>: Correlation vector (C_ij = ⟨s_i s_j⟩ - ⟨s_i⟩⟨s_j⟩)
 */
void saveMagCorr(const string& filename, const vector<double>& si, const vector<double>& sisj, const vector<double>& C) {
    int N_spins = si.size();
    int N_duplet = sisj.size();

    ofstream file(filename.c_str());
    file << N_spins << endl;

    for (int i = 0; i < N_duplet; ++i) {
        file << " " << sisj[i] << " " << C[i];
        if (i < N_spins)
            file << " " << si[i];
        file << endl;
    }

    file.close();
    cout << "File saved: " << filename << endl;
}