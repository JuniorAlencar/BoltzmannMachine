#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <bits/stdc++.h>

using namespace std;

// Function to generate external field values h_i
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

// Function to generate a symmetric matrix J_ij with Gaussian-distributed values
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

// Function to generate a binary spin vector (+1.0 or -1.0)
vector<double> generateBinaryVector(int N) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(0, 1);

    vector<double> vec(N);
    for (int i = 0; i < N; ++i) {
        vec[i] = (dist(gen) == 0) ? -1.0 : 1.0;
    }
    return vec;
}

// Function to compute the Hamiltonian H for a given state
double computeHamiltonian(const vector<double>& sigma, const vector<double>& h, const vector<vector<double>>& J) {
    int N = sigma.size();
    double H = 0.0;

    // External field contribution: sum_i (h_i * sigma_i)
    for (int i = 0; i < N; ++i) {
        H += h[i] * sigma[i];
    }

    // Interaction term: sum_i sum_j>i (J_ij * sigma_i * sigma_j)
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            H -= J[i][j] * sigma[i] * sigma[j];
        }
    }

    return H;
}

// Function to compute the si syntetic
vector<double> computeSi(const vector<vector<double>>& sigma_states) {

    int N_samples = sigma_states.size(); // Number of rows (M_states)
    int N_spins = (N_samples > 0) ? sigma_states[0].size() : 0; // Number of columns (N_spins)
    
    vector<double> s(N_spins, 0.0);
    
    for (int p = 0; p < N_spins; p++) {
        for (int w = 0; w < N_samples; w++) {
            s[p] += sigma_states[w][p]; // Now correctly summing over rows
        }
        s[p] /= N_samples;
    }

    return s;
}

// Function to compute the sisj syntetic
vector<double> computeSiSj(const vector<vector<double>>& sigma_states) {
    int N_samples = sigma_states.size();
    int N_spins = (N_samples > 0) ? sigma_states[0].size() : 0;
    int N_duplet = (N_spins * (N_spins - 1)) / 2;
    
    vector<double> ss(N_duplet, 0.0); // Fix size

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


// Function to save h_i values in a column format
void saveVectorAsColumn(const string& filename, const vector<double>& data) {
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

// Function to save only the upper triangular part of the J matrix
void saveUpperTriangularJMatrix(const string& filename, const vector<vector<double>>& J) {
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

// Function to save multiple sigma states as rows
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

void saveSi(const string& filename, const vector<double>& si){
    ofstream file_si (filename.c_str());
    int N_spins = si.size();
    
    for (int i = 0; i < N_spins; i++)
    {
        file_si << si[i] << endl;
    }

    //fechando arquivos
    file_si.close();
    cout << "File saved: " << filename << endl;
}

void saveSiSj(const string& filename, const vector<double>& sisj){
    ofstream file_sisj (filename.c_str());
    int N_duplet = sisj.size();
    
    for (int i = 0; i < N_duplet; i++)
    {
        file_sisj << sisj[i] << endl;
    }

    //fechando arquivos
    file_sisj.close();
    cout << "File saved: " << filename << endl;
}