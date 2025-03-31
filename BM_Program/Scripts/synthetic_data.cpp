#include "include/synthetic_data.h"

#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <cmath>
#include <numeric>

using namespace std;

#include <cmath>
#include <numeric>

using namespace std;

// Geração de estados com jitter e ruído externo
vector<vector<double>> generateDiverseCorrelationStates(int N, int M, int num_groups, double base_strength, double strength_jitter, double noise_prob) {
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

// Cálculo de si, sisj e C
void computeMagCorr(const vector<vector<double>>& states, vector<double>& si, vector<double>& sisj, vector<double>& C) {
    int M = states.size();
    int N = states[0].size();

    si.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
        for (int m = 0; m < M; ++m)
            si[i] += states[m][i];
        si[i] /= M;
    }

    for (int i = 0; i < N - 1; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double mean_prod = 0.0;
            for (int m = 0; m < M; ++m)
                mean_prod += states[m][i] * states[m][j];
            mean_prod /= M;

            sisj.push_back(mean_prod);
            C.push_back(mean_prod - si[i] * si[j]);
        }
    }
}

// Salvamento em arquivo no formato compatível
void saveMagCorR(const string& filename, const vector<double>& si, const vector<double>& sisj, const vector<double>& C) {
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

int main(int argc, char *argv[]) {

    // Gaussian Parameters
	int n;
	double mean = 1;
	double sigma = 0.25;
	double k = 10;
	int type;
	int H;
    
    string method = argv[1];        // method
    int N_spins = stoi(argv[2]);    // Number of spins
    int M_states = stoi(argv[3]);   // Number of states
    

    vector<double> h = generateHVector(N_spins, -1.0, 1.0);

    // Generate J_ij values
    vector<vector<double>> J = generateJMatrix(N_spins, mean, sigma);
    
    // Generate M sigma states with correlation
    // double target_mean_C = 0.0045;
    // double target_std_C = 0.015;
    // // Generate M sigma states
    // vector<vector<double>> sigmaStates = generateCorrelatedStates(N_spins, M_states, target_mean_C, target_std_C);
    
    int num_groups = 4;
    double base_strength = 0.945;
    double strength_jitter = 0.02;
    double noise_prob = 0.015;
    auto sigmaStates = generateDiverseCorrelationStates(N_spins, M_states, num_groups, base_strength, strength_jitter, noise_prob);
    // Generate M sigma states without correlation
    //vector<vector<double>> sigmaStates = generateUncorrelatedStates(N_spins, M_states);
    
    vector<double> hamiltonianValues = computeHamiltonian(sigmaStates, h, J);

    // Calculate si and sisj and C_i from synthetic data
    vector<double> si_synthetic = computeSi(sigmaStates);
    vector<double> sisj_synthetic = computeSiSj(sigmaStates);
    vector<double> C_synthetic = computeC(si_synthetic, sisj_synthetic);
    
    string file_h_syntetic = "../tests/" + method +  "/hi/h_syntetic.csv";
    string file_si_syntetic = "../tests/" + method +  "/si/si_syntetic.csv";
    string file_sisj_syntetic = "../tests/" + method +  "/sisj/sisj_syntetic.csv";
    string file_j_syntetic = "../tests/" + method +  "/Jij/J_syntetic.csv";
    string file_H_syntetic = "../tests/" + method +  "/H/H_syntetic.csv";
    
    // Save synthetic data
    string file_states_syntetic = "../tests/data_syntetic.dat";
    // Save synthetic mag_corr
    string file_mag_corr_syntetic = "../tests/mag_corr_syntetic.dat";
    
    // Save h_i values as a column
    saveVectorAsColumn(file_h_syntetic, h);

    // Save upper triangular J_ij matrix
    saveUpperTriangularJMatrix(file_j_syntetic, J);

    // Save sigma states (each row = one state)
    saveSigmaStates(file_states_syntetic, sigmaStates);

    // Save Hamiltonian values (each row corresponds to a sigma state)
    saveVectorAsColumn(file_H_syntetic, hamiltonianValues);

    saveSi(file_si_syntetic, si_synthetic);
    saveSiSj(file_sisj_syntetic, sisj_synthetic);
    //saveMagCorr(file_mag_corr_syntetic, si_synthetic, sisj_synthetic, C_synthetic);
    

    computeMagCorr(sigmaStates, si_synthetic, sisj_synthetic, C_synthetic);
    saveMagCorR("../tests/mag_corr_fully_negative_aligned.dat", si_synthetic, sisj_synthetic, C_synthetic);

    // Calcular média e desvio padrão de C
    double mean_C = 0.0, std_C = 0.0;
    for (double c : C_synthetic) mean_C += c;
    mean_C /= C_synthetic.size();
    
    for (double c : C_synthetic) std_C += (c - mean_C) * (c - mean_C);
    std_C = sqrt(std_C / C_synthetic.size());
    
    cout << "Correlação média: " << mean_C << " ± " << std_C << endl;
    
    // Saving information about synthetic data
    string info = "../tests/" + method + "/info.dat";
    ofstream file_info (info.c_str());
    file_info << "N_spins = " << N_spins << "\tN_samples = " << M_states << endl;
    file_info.close();

    return 0;
}
