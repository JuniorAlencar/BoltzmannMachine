#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <bits/stdc++.h>
#include <cmath>
#include <numeric>
#include <optional>
#include <queue>
#include "nr3.h"
#include "network.h"

using namespace std;
using std::mt19937;

// If you use vscode, install the Doxygen Documentation Generator extension

// Compute properties =========================================================/

// /**
//  * @brief Generates a 1D vector containing Gaussian-distributed J_ij values.
//  * 
//  * @param N_pairs Number of spin pairs (should be N*(N-1)/2 for N spins).
//  * @param stddev Standard deviation of the normal distribution.
//  * @param gen Reference to a random number generator.
//  * @param mean Mean of the normal distribution (default = 0.0).
//  * @return vector<double> Vector containing J_ij values.
//  */
// vector<double> ComputeJValues(int N_pairs, double stddev, mt19937 &gen, double mean = 0.0) {
//     normal_distribution<double> dist(mean, stddev);

//     vector<double> J(N_pairs, 0.0);
//     for (int i = 0; i < N_pairs; i++)
//         J[i] = dist(gen);

//     return J;
// }


// /**
//  * @brief Generates a vector of external field values h_i sampled from a uniform distribution.
//  * 
//  * @param N Number of spins.
//  * @param min Minimum value of the uniform distribution (e.g., -1.0).
//  * @param max Maximum value of the uniform distribution (e.g., 1.0).
//  * @param gen Reference to a random number generator (mt19937).
//  * @return vector<double> Vector containing the h_i values.
//  */
// vector<double> ComputehValues(int N, double min, double max, mt19937 &gen) {
//     uniform_real_distribution<double> dist(min, max);

//     vector<double> h(N);
//     for (int i = 0; i < N; ++i) {
//         h[i] = dist(gen);
//     }

//     return h;
// }


#include <vector>
#include <random>
#include <cmath>
#include <queue>

using namespace std;

void wolff_update(Rede &rede, double beta, mt19937 &gen) {
    int N = rede.n;
    int N_pairs = (N * (N - 1)) / 2;

    uniform_real_distribution<double> dist(0.0, 1.0);
    uniform_int_distribution<int> spin_dist(0, N - 1);

    // Probabilidade de adicionar um vizinho ao cluster
    auto P_add = [&](double Jij) {
        return 1.0 - exp(-2.0 * beta * fabs(Jij));
    };

    // Inicializar o cluster
    vector<bool> in_cluster(N, false);

    // Escolher spin raiz aleatório
    int root_spin = spin_dist(gen);
    int s_value = rede.s[root_spin];
    in_cluster[root_spin] = true;

    queue<int> cluster_queue;
    cluster_queue.push(root_spin);

    // Crescer o cluster
    while (!cluster_queue.empty()) {
        int i = cluster_queue.front();
        cluster_queue.pop();

        int idx = 0;
        for (int j = 0; j < N - 1; ++j) {
            for (int k = j + 1; k < N; ++k, ++idx) {
                int a = j, b = k;
                if (a == i || b == i) {
                    int neighbor = (a == i) ? b : a;

                    if (!in_cluster[neighbor] && rede.s[neighbor] == s_value) {
                        double Jij = rede.J[idx];
                        if (dist(gen) < P_add(Jij)) {
                            in_cluster[neighbor] = true;
                            cluster_queue.push(neighbor);
                        }
                    }
                }
            }
        }
    }

    // Flipar o cluster inteiro
    for (int i = 0; i < N; ++i) {
        if (in_cluster[i]) {
            rede.s[i] *= -1;
        }
    }
}

void GenerateStates_Wolff(
    Rede &rede,
    vector<double> &av_s,
    vector<double> &av_ss,
    const int t_eq,
    const int relx,
    const int rept,
    const int M,
    const double beta,
    vector<double> &energies,
    vector<vector<int>> &states,
    mt19937 &gen
) {
    int N = rede.n;
    int N_pairs = (N * (N - 1)) / 2;

    av_s.assign(N, 0.0);
    av_ss.assign(N_pairs, 0.0);
    energies.clear();
    states.assign(M, vector<int>(N, 0));

    // ---------------- Equilibration ----------------
    for (int t = 0; t < t_eq; ++t)
        wolff_update(rede, beta, gen);

    int m_count = 0;

    for (int rep = 0; rep < rept; ++rep) {
        while (m_count < M) {

            // Para cada amostra:
            int relax_count = 0;
            while (relax_count < relx) {
                wolff_update(rede, beta, gen);
                relax_count++;
            }

            // Coleta de estado
            // (skip automático: se estado atual é igual ou invertido demais do anterior, pula)
            bool accept = true;

            if (m_count > 0) {
                double corr = 0.0;
                for (int k = 0; k < N; ++k)
                    corr += states[m_count - 1][k] * rede.s[k];

                corr /= N;

                // Se a correlação for muito próxima de -1 (flip completo), não coleta
                if (fabs(corr + 1.0) < 1e-3)  // tolerância ajustável
                    accept = false;
            }

            if (!accept) continue; // Pula e faz próxima relaxação

            // ---- Calcular energia ----
            double energy = 0.0;
            int idx = 0;
            for (int i = 0; i < N; ++i)
                energy -= rede.h[i] * rede.s[i];
            for (int i = 0; i < N - 1; ++i)
                for (int j = i + 1; j < N; ++j)
                    energy -= rede.J[idx++] * rede.s[i] * rede.s[j];

            energies.push_back(energy);

            // ---- Acumular observáveis ----
            for (int k = 0; k < N; ++k)
                av_s[k] += rede.s[k];

            idx = 0;
            for (int i = 0; i < N - 1; ++i)
                for (int j = i + 1; j < N; ++j)
                    av_ss[idx++] += rede.s[i] * rede.s[j];

            // ---- Salvar estado ----
            for (int k = 0; k < N; ++k)
                states[m_count][k] = rede.s[k];

            ++m_count;
        }
    }

    // ---- Normalizar médias ----
    for (int i = 0; i < N; ++i)
        av_s[i] /= M;

    for (int i = 0; i < N_pairs; ++i)
        av_ss[i] /= M;
}



/**
 * @brief Gera estados e momentos <s_i> e <s_i s_j> para uma rede de Ising usando Metropolis.
 *
 * @param rede      Objeto Rede com spins, J, h
 * @param av_s      Vetor para armazenar as médias <s_i>
 * @param av_ss     Vetor para armazenar as médias <s_i s_j>
 * @param t_eq      Passos de equilíbrio (burn-in)
 * @param relx      Relaxamento entre amostras
 * @param rept      Número de repetições do processo
 * @param M         Número total de amostras desejadas
 * @param beta      Inverso da temperatura
 * @param energies  Vetor para armazenar energias
 * @param states    Vetor para armazenar estados de spin
 * @param gen       Semente randômica
 */
void GenerateStates(
    Rede &rede,
    vector<double> &av_s,
    vector<double> &av_ss,
    const int t_eq,
    const int relx,
    const int rept,
    const int M,
    const double beta,
    vector<double> &energies,
    vector<vector<int>> &states,
    mt19937 &gen
) {
    int N = rede.n;
    int N_pairs = (N * (N - 1)) / 2;
    
    uniform_real_distribution<double> dist(0.0, 1.0);
    uniform_int_distribution<int> spin_dist(0, N - 1);

    av_s.assign(N, 0.0);
    av_ss.assign(N_pairs, 0.0);
    energies.clear();
    states.assign(M, vector<int>(N, 0));

    int m_count = 0;

    for (int rep = 0; rep < rept; ++rep) {
        while (m_count < M) {
            // Equilibration
            for (int t = 0; t < t_eq; ++t) {
                int i = spin_dist(gen);
                double dE = 2.0 * rede.s[i] * rede.h[i];
                int idx = 0;
                for (int a = 0; a < N - 1; ++a)
                    for (int b = a + 1; b < N; ++b, ++idx)
                        if (a == i)
                            dE += 2.0 * rede.s[i] * rede.s[b] * rede.J[idx];
                        else if (b == i)
                            dE += 2.0 * rede.s[i] * rede.s[a] * rede.J[idx];

                if (dE <= 0.0 || dist(gen) < exp(-beta * dE)) {
                    rede.s[i] *= -1;

                    // Save energy after accepted flip
                    double energy = 0.0;
                    int idx3 = 0;
                    for (int a = 0; a < N; ++a)
                        energy -= rede.h[a] * rede.s[a];
                    for (int a = 0; a < N - 1; ++a)
                        for (int b = a + 1; b < N; ++b)
                            energy -= rede.J[idx3++] * rede.s[a] * rede.s[b];
                    energies.push_back(energy);
                }
            }

            // Sampling
            for (int relax = 0; relax < relx && m_count < M; ++relax) {
                int i = spin_dist(gen);
                double dE = 2.0 * rede.s[i] * rede.h[i];
                int idx = 0;
                for (int a = 0; a < N - 1; ++a)
                    for (int b = a + 1; b < N; ++b, ++idx)
                        if (a == i)
                            dE += 2.0 * rede.s[i] * rede.s[b] * rede.J[idx];
                        else if (b == i)
                            dE += 2.0 * rede.s[i] * rede.s[a] * rede.J[idx];

                if (dE <= 0.0 || dist(gen) < exp(-beta * dE)) {
                    rede.s[i] *= -1;

                    // Save energy after accepted flip
                    double energy = 0.0;
                    int idx3 = 0;
                    for (int a = 0; a < N; ++a)
                        energy -= rede.h[a] * rede.s[a];
                    for (int a = 0; a < N - 1; ++a)
                        for (int b = a + 1; b < N; ++b)
                            energy -= rede.J[idx3++] * rede.s[a] * rede.s[b];
                    energies.push_back(energy);
                }
            }

            // Accumulate observables
            for (int k = 0; k < N; ++k)
                av_s[k] += rede.s[k];

            int idx2 = 0;
            for (int a = 0; a < N - 1; ++a)
                for (int b = a + 1; b < N; ++b)
                    av_ss[idx2++] += rede.s[a] * rede.s[b];

            // Save spin configuration
            for (int k = 0; k < N; ++k)
                states[m_count][k] = rede.s[k];

            ++m_count;
        }
    }

    // Normalize moments
    for (int i = 0; i < N; ++i)
        av_s[i] /= M;
    for (int i = 0; i < N_pairs; ++i)
        av_ss[i] /= M;
}


/**
 * @brief Computes the first moment (magnetization) of the system.
 *
 * @param sigma_states \c vector<vector<double>>: Matrix containing all spin states (M samples × N spins)
 * @return \c vector<double>: Magnetization vector s_i
 */
vector<double> computeSi(const vector<vector<int>>& sigma_states) {
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
vector<double> computeSiSj(const vector<vector<int>>& sigma_states) {
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
 * @brief Computes the third-order moment ⟨s_i s_j s_k⟩ for all i < j < k.
 *
 * @param sigma_states Vector of spin configurations (M × N).
 * @return Flattened vector of size C(N, 3) with ⟨s_i s_j s_k⟩ values.
 */
std::vector<double> computeSiSjSk(const std::vector<std::vector<int>>& sigma_states) {
    int M = sigma_states.size();
    int N = (M > 0) ? sigma_states[0].size() : 0;
    std::vector<double> sss;

    for (int i = 0; i < N - 2; ++i) {
        for (int j = i + 1; j < N - 1; ++j) {
            for (int k = j + 1; k < N; ++k) {
                double acc = 0.0;
                for (int m = 0; m < M; ++m)
                    acc += sigma_states[m][i] * sigma_states[m][j] * sigma_states[m][k];
                sss.push_back(acc / M);
            }
        }
    }
    return sss;
}

/**
 * @brief Computes Pearson correlation coefficients P_{ij} = C_{ij} / (σ_i σ_j) in flat format.
 *
 * @param sigma_states Vector of spin configurations (M × N).
 * @param cov Vector of covariance values from computeCovariance.
 * @return Flattened vector of Pearson correlation values P_{ij} for i < j.
 */
std::vector<double> computePij(const std::vector<std::vector<int>>& sigma_states, const std::vector<double>& cov) {
    int M = sigma_states.size();
    int N = (M > 0) ? sigma_states[0].size() : 0;

    std::vector<double> mean(N, 0.0);
    std::vector<double> stddev(N, 0.0);
    std::vector<double> Pij(cov.size(), 0.0);

    for (int i = 0; i < N; ++i)
        for (int m = 0; m < M; ++m)
            mean[i] += sigma_states[m][i];
    for (int i = 0; i < N; ++i)
        mean[i] /= M;

    for (int i = 0; i < N; ++i)
        for (int m = 0; m < M; ++m)
            stddev[i] += std::pow(sigma_states[m][i] - mean[i], 2);
    for (int i = 0; i < N; ++i)
        stddev[i] = std::sqrt(stddev[i] / M);

    int idx = 0;
    for (int i = 0; i < N - 1; ++i) {
        for (int j = i + 1; j < N; ++j, ++idx) {
            double denom = stddev[i] * stddev[j];
            Pij[idx] = (denom != 0.0) ? cov[idx] / denom : 0.0;
        }
    }

    return Pij;
}

/**
 * @brief Computes third-order central moments T_{ijk} = ⟨s_i s_j s_k⟩ - s_i⟨s_j s_k⟩ - ... + 2s_i s_j s_k.
 *
 * @param sigma_states Vector of spin configurations (M × N).
 * @param third_moment Precomputed third-order moments from computeThirdMoment.
 * @return Flattened vector of central moment values T_{ijk}.
 */
std::vector<double> computeTriplet(const std::vector<std::vector<int>>& sigma_states, const std::vector<double>& third_moment) {
    int M = sigma_states.size();
    int N = (M > 0) ? sigma_states[0].size() : 0;
    std::vector<double> s(N, 0.0);
    std::vector<double> Tijk;

    for (int i = 0; i < N; ++i)
        for (int m = 0; m < M; ++m)
            s[i] += sigma_states[m][i];
    for (int i = 0; i < N; ++i)
        s[i] /= M;

    int idx = 0;
    for (int i = 0; i < N - 2; ++i) {
        for (int j = i + 1; j < N - 1; ++j) {
            for (int k = j + 1; k < N; ++k) {
                double sij = 0.0, sik = 0.0, sjk = 0.0;
                for (int m = 0; m < M; ++m) {
                    sij += sigma_states[m][i] * sigma_states[m][j];
                    sik += sigma_states[m][i] * sigma_states[m][k];
                    sjk += sigma_states[m][j] * sigma_states[m][k];
                }
                sij /= M; sik /= M; sjk /= M;

                double T = third_moment[idx] - s[i]*sjk - s[j]*sik - s[k]*sij + 2*s[i]*s[j]*s[k];
                Tijk.push_back(T);
                idx++;
            }
        }
    }

    return Tijk;
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
            C[aux++] = ss[aux] - s[p] * s[pp];
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
vector<double> computeHamiltonian(const vector<vector<int>>& sigmaStates, const vector<double>& h, const vector<double>& J) {
    int M = sigmaStates.size();
    int N = (M > 0) ? sigmaStates[0].size() : 0;
    vector<double> H_values(M, 0.0);
    
    for (int m = 0; m < M; ++m) {
        double H = 0.0;

        // External field contribution
        for (int i = 0; i < N; ++i) {
            H += h[i] * sigmaStates[m][i];
        }
        int index = 0;
        // Pairwise interaction term (only upper triangle considered)
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                H -= J[index] * sigmaStates[m][i] * sigmaStates[m][j];
            }
            index ++;
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
void computeMagCorr(const vector<vector<int>>& sigmaStates, vector<double>& si, vector<double>& sisj, vector<double>& C) {
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
 * @param filename \c string: Path + filename (including extension)
 * @param data \c vector<double>: Vector containing h_i values
 */
void saveValues(const string& filename,const string& header ,const vector<double>& data) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error opening file " << filename << endl;
        return;
    }
    file << header << endl;
    for (const auto& value : data) {
        file << fixed << setprecision(10) << value << "\n";
    }
    file.close();
    cout << "File saved: " << filename << endl;
}


/**
 * @brief Saves multiple spin configurations (sigma states) row-wise in CSV-like format.
 *
 * @param filename \c string: Path + filename (including extension)
 * @param sigmaStates \c vector<vector<double>>: Matrix with spin states (each row = one configuration)
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

void saveSigmaStates(const string& filename, const vector<vector<int>>& sigmaStates) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error opening file " << filename << endl;
        return;
    }
    int N = sigmaStates[0].size(); // Number of spins

    // Write CSV header
    for (int i = 1; i <= N; ++i) {
        file << "s_" << i;
        if (i < N)
            file << ",";
    }
    file << "\n";

    // Write data with one decimal place, no explicit positive sign
    for (const auto& state : sigmaStates) {
        for (size_t i = 0; i < state.size(); ++i) {
            file << fixed << setprecision(1) << static_cast<double>(state[i]);
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
 * @param filename \c string: Name of file, including path and extension
 * @param si \c vector<double>: Vector containing the mean value of each spin (s_i)
 */
void saveS(const string& filename, const string& header ,const vector<double>& si) {
    ofstream file_si(filename.c_str());
    int N_spins = si.size();
    file_si << header << endl;
    for (int i = 0; i < N_spins; i++) {
        file_si << si[i] << endl;
    }

    file_si.close();
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

/**
 * @brief Computes the total energy of a spin configuration.
 * 
 * @param state Vector of spin states (+1 or -1) for each site.
 * @param h Vector of external field values h_i.
 * @param J Flattened upper-triangular vector of interaction values J_ij (i < j).
 * @return double Energy of the system.
 */
void SaveEnergy(const string& filename, const vector<double>& energies) {
    int N_steps = energies.size();
    
    ofstream file(filename.c_str());
    
    file << "time" << " " << "energy_state" << endl;
    for(int i=0; i < N_steps; i++){
        file << i << " " << energies[i] << endl;
    }
    file.close();
    cout << "File saved: " << filename << endl;
}