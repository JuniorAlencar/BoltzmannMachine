#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <bits/stdc++.h>

using namespace std;


// Compute properties =========================================================/

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

vector<double> computeC(const vector<double> s, const vector<double> ss){
    int N_spins = s.size();
    int N_duplet = ss.size();
    
    vector<double> C(N_duplet, 0.0);
    int aux = 0;
    for (int p = 0; p < N_spins - 1; p++)
	{
		for (int pp = p+1; pp < N_spins; pp++)
		{
			C[aux] = ss[aux] - s[p] * s[pp];			
			aux++;
		}
	}
    return C;
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

vector<vector<double>> generateCorrelatedStates(int N, int M,const double target_mean_C, const double target_std_C) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    normal_distribution<double> corrDist(target_mean_C, target_std_C);
    
    vector<vector<double>> sigmaStates(M, vector<double>(N));
    
    // Gerar primeiro estado aleatório
    for (int i = 0; i < N; ++i) {
        sigmaStates[0][i] = (dist(gen) < 0.5) ? -1.0 : 1.0;
    }
    
    // Criar uma matriz de correlação alvo
    vector<double> target_C((N * (N - 1)) / 2);
    for (double& c : target_C) {
        c = corrDist(gen);
    }
    
    // Gerar estados subsequentes tentando obedecer a correlação alvo
    int index = 0;
    for (int m = 1; m < M; ++m) {
        for (int i = 0; i < N; ++i) {
            sigmaStates[m][i] = (dist(gen) < 0.5) ? -1.0 : 1.0;
        }
        
        // Ajuste para manter a correlação esperada
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

vector<vector<double>> generateUncorrelatedStates(int N, int M) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    
    vector<vector<double>> sigmaStates(M, vector<double>(N));
    
    for (int m = 0; m < M; ++m) {
        for (int i = 0; i < N; ++i) {
            sigmaStates[m][i] = (dist(gen) < 0.5) ? -1.0 : 1.0;
        }
    }

    return sigmaStates;
}

// Função para calcular a Hamiltoniana para cada estado
vector<double> computeHamiltonian(const vector<vector<double>>& sigmaStates, const vector<double>& h, const vector<vector<double>>& J) {
    int M = sigmaStates.size();
    int N = (M > 0) ? sigmaStates[0].size() : 0;
    vector<double> H_values(M, 0.0);
    
    for (int m = 0; m < M; ++m) {
        double H = 0.0;
        
        // Contribuição do campo externo
        for (int i = 0; i < N; ++i) {
            H += h[i] * sigmaStates[m][i];
        }
        
        // Termo de interação
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                H -= J[i][j] * sigmaStates[m][i] * sigmaStates[m][j];
            }
        }
        
        H_values[m] = H;
    }
    return H_values;
}


// Save properties =========================================================/

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

void saveMagCorr(const string& filename, const vector<double> si, const vector<double>& sisj, const vector<double> C){
    int N_spins = si.size();
    int N_duplet = sisj.size();
    
    ofstream file_MagCorr (filename.c_str());
    
    
    file_MagCorr << N_spins << endl;
	
	for (int i = 0; i < N_duplet; i++)
	{
		file_MagCorr << " " << sisj[i] << " " << C[i];
		
		if (i < N_spins)
			file_MagCorr << " " << si[i];
			
		file_MagCorr << endl;
	}

    file_MagCorr.close();
    cout << "File saved: " << filename << endl;
}