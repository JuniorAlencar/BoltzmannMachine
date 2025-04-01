#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <bits/stdc++.h>
#include <cmath>
#include <numeric>

using namespace std;

// If you use vscode, install the Doxygen Documentation Generator extension

// Compute properties =========================================================/

/**
 * @brief Generates a symmetric matrix J_ij with Gaussian-distributed values.
 *
 * @param N **int**: Number of spins
 * @param mean **double**: Mean of the Gaussian distribution
 * @param stddev **double**: Standard deviation of the Gaussian distribution
 * @return Matrix J of size NxN with symmetric Gaussian-distributed values
 */
vector<vector<double>> generateJMatrix(int N, double mean, double stddev) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist(mean, stddev);

    vector<vector<double>> J(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) { // Fill only upper triangle
            J[i][j] = J[j][i] = dist(gen);
        }
    }
    return J;
}

/**
 * @brief Generates M correlated spin states with approximately desired pairwise correlations.
 *
 * @param N **int**: Number of spins
 * @param M **int**: Number of states (samples)
 * @param target_mean_C **double**: Desired average pairwise correlation
 * @param target_std_C **double**: Desired standard deviation of pairwise correlations
 * @return Matrix of size MxN where each row is a spin configuration
 */
vector<vector<double>> generateCorrelatedStates(int N, int M, const double target_mean_C, const double target_std_C) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    normal_distribution<double> corrDist(target_mean_C, target_std_C);
    
    vector<vector<double>> sigmaStates(M, vector<double>(N));
    
    // Generate first random state
    for (int i = 0; i < N; ++i) {
        sigmaStates[0][i] = (dist(gen) < 0.5) ? -1.0 : 1.0;
    }
    
    // Target correlation values
    vector<double> target_C((N * (N - 1)) / 2);
    for (double& c : target_C) {
        c = corrDist(gen);
    }
    
    // Generate subsequent states trying to match target correlations
    int index = 0;
    for (int m = 1; m < M; ++m) {
        for (int i = 0; i < N; ++i) {
            sigmaStates[m][i] = (dist(gen) < 0.5) ? -1.0 : 1.0;
        }

        index = 0;
        for (int i = 0; i < N - 1; i++) {
            for (int j = i + 1; j < N; j++) {
                if ((dist(gen) < (target_C[index] + 1.0) / 2.0)) {
                    sigmaStates[m][j] = sigmaStates[m][i];
                }
                index++;
            }
        }
    }
    
    return sigmaStates;
}

/**
 * @brief Generates a vector of external field values h_i from a uniform distribution.
 *
 * @param N **int**: Number of spins
 * @param min **double**: Minimum value of the external field
 * @param max **double**: Maximum value of the external field
 * @return Vector of size N with uniformly distributed values in [min, max]
 */
vector<double> generateHVector(int N, double min, double max) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(min, max);

    vector<double> vec(N);
    for (int i = 0; i < N; ++i) {
        vec[i] = dist(gen);
    }
    return vec;
}



/**
 * @brief Function to generate states simulating real physical systems (without Monte Carlo)
 * 
 * @param N Number of spins **(int)**
 * @param M Number of states (samples) **(int)**
 * @param num_groups Number of groups into which the N spins are divided **(int)**
 * @param base_strength Base probability of each spin following the leader of its group **(double)**
 * @param strength_jitter Random variation in "base strength" for each group **(double)**
 * @param noise_prob Probability of each spin inverting independently of its group **(double)**
* Example:
* ```
* num_groups = 4; 
* base_strength = 0.945;
* strength_jitter = 0.01;
* noise_prob = 0.005;
 * ```
 * @return Matrix of spin states (vector<vector<double>>)
 **/
vector<vector<double>> generateDiverseCorrelationStates(int N, int M, int num_groups, 
                                                        double base_strength, double strength_jitter, 
                                                        double noise_prob) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);

    vector<vector<double>> sigmaStates(M, vector<double>(N));

    // Criar grupos
    vector<vector<int>> group_indices(num_groups);
    int spins_per_group = N / num_groups;
    for (int g = 0; g < num_groups; ++g)
        for (int i = g * spins_per_group; i < (g + 1) * spins_per_group; ++i)
            group_indices[g].push_back(i);

    for (int i = spins_per_group * num_groups; i < N; ++i)
        group_indices.back().push_back(i);

    // Força de grupo com jitter
    vector<double> group_strengths(num_groups);
    for (int g = 0; g < num_groups; ++g)
        group_strengths[g] = max(0.0, min(1.0, base_strength + (dist(gen) * 2 - 1) * strength_jitter));

    // Geração dos estados
    for (int m = 0; m < M; ++m) {
        for (int g = 0; g < num_groups; ++g) {
            for (int idx : group_indices[g]) {
                bool aligned = dist(gen) < group_strengths[g];
                sigmaStates[m][idx] = aligned ? -1.0 : 1.0;
            }
        }

        // Ruído aleatório
        for (int i = 0; i < N; ++i) {
            if (dist(gen) < noise_prob)
                sigmaStates[m][i] *= -1.0;
        }
    }

    return sigmaStates;
}

/**
 * @brief Computes the first moment (magnetization) of the system.
 *
 * @param sigma_states vector<vector<double>>: Matrix containing all spin states (M samples × N spins)
 * @return vector<double>: Magnetization vector s_i
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
 * @param sigma_states vector<vector<double>>: Matrix containing all spin states (M samples × N spins)
 * @return vector<double>: Vector of s_i * s_j averages for all unique pairs (i < j)
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
 * @param s vector<double>: Magnetization vector ⟨s_i⟩
 * @param ss vector<double>: Second moment vector ⟨s_i s_j⟩
 * @return vector<double>: Connected correlation values C_ij
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
 * @param sigmaStates vector<vector<double>>: Matrix of spin states (M samples × N spins)
 * @param h vector<double>: External field vector h_i
 * @param J vector<vector<double>>: Symmetric interaction matrix J_ij
 * @return vector<double>: Vector of Hamiltonian values for each configuration
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
 * @param sigmaStates vector<vector<double>>: Matrix containing all the spin states (M samples × N spins)
 * @param si vector<double>&: Output vector of mean spin values (first moment)
 * @param sisj vector<double>&: Output vector of mean products s_i * s_j (second moment)
 * @param C vector<double>&: Output vector of connected correlations C_ij
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
 * @param filename string: Path + filename (including extension)
 * @param data vector<double>: Vector containing h_i values
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
 * @param filename string: Path + filename (including extension)
 * @param J vector<vector<double>>: Symmetric interaction matrix J_ij
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
 * @param filename string: Path + filename (including extension)
 * @param sigmaStates vector<vector<double>>: Matrix with spin states (each row = one configuration)
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
 * @param filename string: Name of file, including path and extension
 * @param si vector<double>: Vector containing the mean value of each spin (s_i)
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
 * @param filename string: Name of file, including path and extension
 * @param sisj vector<double>: Vector of pairwise spin products
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
 * Example:
 * ```
 * 4
 * 0.82 0.05 0.91
 * 0.76 0.03
 * 0.88 0.10 0.95
 * ```
 *
 * @param filename string: Name of file, including path and extension
 * @param si vector<double>: First moment (mean values of spins)
 * @param sisj vector<double>: Second moment (⟨s_i s_j⟩)
 * @param C vector<double>: Correlation vector (C_ij = ⟨s_i s_j⟩ - ⟨s_i⟩⟨s_j⟩)
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