#include "include/synthetic_data.h"

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
    
    // Generate h_i values
    vector<double> h = generateHVector(N_spins, -1.0, 1.0);

    // Generate J_ij values
    vector<vector<double>> J = generateJMatrix(N_spins, mean, sigma);

    double target_mean_C = 0.0045;
    double target_std_C = 0.015;
    // Generate M sigma states
    
    vector<vector<double>> sigmaStates = generateCorrelatedStates(N_spins, M_states, target_mean_C, target_std_C);
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
    saveMagCorr(file_mag_corr_syntetic, si_synthetic, sisj_synthetic, C_synthetic);

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
