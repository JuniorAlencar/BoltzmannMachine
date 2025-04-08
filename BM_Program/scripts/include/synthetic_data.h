#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <bits/stdc++.h>
#include <cmath>
#include <numeric>
#include <optional>

using namespace std;

// If you use vscode, install the Doxygen Documentation Generator extension

// Compute properties =========================================================/


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
        int step = 0;
        while (m_count < M) {
            // Equilíbrio + relaxamento
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
                if (dE <= 0.0 || dist(gen) < exp(-beta * dE))
                    rede.s[i] *= -1;
            }

            // Coleta
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
                if (dE <= 0.0 || dist(gen) < exp(-beta * dE))
                    rede.s[i] *= -1;
            }

            // Acumular momentos
            for (int k = 0; k < N; ++k)
                av_s[k] += rede.s[k];

            int idx2 = 0;
            for (int a = 0; a < N - 1; ++a)
                for (int b = a + 1; b < N; ++b)
                    av_ss[idx2++] += rede.s[a] * rede.s[b];

            // Energia
            double energy = 0.0;
            int idx3 = 0;
            for (int a = 0; a < N; ++a)
                energy -= rede.h[a] * rede.s[a];
            for (int a = 0; a < N - 1; ++a)
                for (int b = a + 1; b < N; ++b)
                    energy -= rede.J[idx3++] * rede.s[a] * rede.s[b];
            energies.push_back(energy);
            cout << "update state" + to_string(m_count) << endl;
            for (int k = 0; k < N; ++k)
                states[m_count][k] = rede.s[k];

            ++m_count;
        }
    }

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
void saveValues(const string& filename, const vector<double>& data) {
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
 * @brief Saves multiple spin configurations (sigma states) row-wise in CSV-like format.
 *
 * @param filename \c string: Path + filename (including extension)
 * @param sigmaStates \c vector<vector<double>>: Matrix with spin states (each row = one configuration)
 */
void saveSigmaStates(const string& filename, const vector<vector<int>>& sigmaStates) {
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
 * @param filename \c string: Name of file, including path and extension
 * @param si \c vector<double>: Vector containing the mean value of each spin (s_i)
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
 * @param filename \c string: Name of file, including path and extension
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