#include "include/synthetic_data.h"
#include <algorithm>

int main(int argc, char *argv[]) {

    // Gaussian Parameters
	int n;
	double mean = 1;
	double sigma = 0.25;
	double k = 10;
	int type;
	int H;
    
    // Lista de métodos válidos
    vector<string> valid_methods = {
        "exact", "metropolis", "parallel_tempering", "swendsen_wang", "wang_landau"
    };

    if (argc < 4) {
        cerr << "Uso: " << argv[0] << " <method> <N_spins> <M_states>" << endl;
        cerr << "Métodos disponíveis: exact, metropolis, parallel_tempering, swendsen_wang, wang_landau" << endl;
        return 1;
    }

    string method = argv[1];

    // Verifica se o método fornecido é válido
    if (find(valid_methods.begin(), valid_methods.end(), method) == valid_methods.end()) {
        cerr << "Erro: método '" << method << "' inválido." << endl;
        cerr << "Métodos válidos: exact, metropolis, parallel_tempering, swendsen_wang, wang_landau" << endl;
        return 1;
    }

    int N_spins = stoi(argv[2]);
    int M_states = stoi(argv[3]);

    cout << "Método selecionado: " << method << endl;
    cout << "Número de spins: " << N_spins << ", Número de estados: " << M_states << endl;
    

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
    double strength_jitter = 0.01;
    double noise_prob = 0.005;
    auto sigmaStates = generateDiverseCorrelationStates(N_spins, M_states, num_groups, base_strength, strength_jitter, noise_prob);
    
    // Generate M sigma states without correlation
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
    string file_states_syntetic = "../tests/data_synteticN" + to_string(N_spins) + ".dat";
    // Save synthetic mag_corr
    string file_mag_corr_syntetic = "../tests/mag_corr_synteticN" + to_string(N_spins) +  ".dat";
    
    // Save h_i values as a column
    savehH(file_h_syntetic, h);

    // Save upper triangular J_ij matrix
    saveJ(file_j_syntetic, J);

    // Save sigma states (each row = one state)
    saveSigmaStates(file_states_syntetic, sigmaStates);

    // Save Hamiltonian values (each row corresponds to a sigma state)
    savehH(file_H_syntetic, hamiltonianValues);

    saveSi(file_si_syntetic, si_synthetic);
    saveSiSj(file_sisj_syntetic, sisj_synthetic);
    //saveMagCorr(file_mag_corr_syntetic, si_synthetic, sisj_synthetic, C_synthetic);
    
    computeMagCorr(sigmaStates, si_synthetic, sisj_synthetic, C_synthetic);
    saveMagCorr(file_mag_corr_syntetic, si_synthetic, sisj_synthetic, C_synthetic);

    // Saving information about synthetic data
    string info = "../tests/" + method + "/info.dat";
    ofstream file_info (info.c_str());
    file_info << "N_spins = " << N_spins << "\tN_samples = " << M_states << endl;
    file_info.close();

    return 0;
}
