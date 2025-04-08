//#include "./include/nr3.h"
#include "./include/network_jr.h"
#include "./include/synthetic_data.h"

#include <algorithm>

int main(int argc, char *argv[]) {

    // Gaussian Parameters
	int n;
	double mean = 0.0;
	double sigma = 0.25;
	double k = 10;
	int type;
	
    int N_spins = stoi(argv[1]);
    int M_states = stoi(argv[2]);
    int seed = stoi(argv[3]);

    if (argc < 4) {
        cerr << "Uso: " << argv[0] << "<N_spins> <M_states> <seed>" << endl;
        
        return 1;
    }
    cout << "Number of spins: " << N_spins << ", Number of states: " << M_states << ", seed: " << seed << endl;
    
    // Generate to code
    std::mt19937 gen(seed);
    
    // Number of combinations in pairs s_i * s_j (s_i * s_j = s_j * s_i) -> symmetric matrix
    int N_pairs = (N_spins * (N_spins - 1)) / 2;
    
    // Initial network, with J and h obtaining from gaussian and uniform distributions, respectively, with mean=0.0
    Rede net(N_spins, 0.0, 1.0, k, type, 1, seed);
    double min = -1.0, max = 1.0;
    vector<double> J_ij = ComputeJValues(N_pairs, sigma, gen, mean);
    vector<double> h_i = ComputehValues(N_spins, min, max, gen);
    net.h = h_i;
    net.J = J_ij;
    
    vector<double> av_s(N_spins, 0.0), av_ss(N_pairs, 0.0);
    
    // Vector to alocate the energies throughout the implementation
    vector<double> energies;
    vector<vector<int>> sigmaStates;
    
    // MONTE CARLO UPDATE ----------------------
    // Monte Carlo variables
	int t_eq = 500; // 150
	int relx = 2500;
	int rept = 40;
    
    // Number of steps
	int t_step = (M_states * relx) / rept;
    
    // inter_max - inter = Number of interations
    int inter = 1;
	int inter_max = 300000;
	
    GenerateStates(net, av_s, av_ss, t_eq,  relx, rept, M_states, 1.0, energies, sigmaStates, gen);
    
    // Compute Hamiltonian
    vector<double> hamiltonianValues = computeHamiltonian(sigmaStates, net.h, net.J);

    // Calculate si and sisj and C_i from synthetic data
    vector<double> si_synthetic = computeSi(sigmaStates);
    vector<double> sisj_synthetic = computeSiSj(sigmaStates);
    vector<double> Cij_synthetic = computeC(si_synthetic, sisj_synthetic);
    vector<double> Pij_synthetic = computePij(sigmaStates, Cij_synthetic);
    vector<double> sisjsk_synthetic = computeSiSjSk(sigmaStates);
    vector<double> Tijk_synthetic = computeTriplet(sigmaStates, sisjsk_synthetic);

    string file_h_syntetic = "../tests/synthetic/hi/h_synteticN" + to_string(N_spins) + ".dat";
    string file_H_syntetic = "../tests/synthetic/H/H_synteticN" + to_string(N_spins) + ".dat";
    string file_j_syntetic = "../tests/synthetic/Jij/J_synteticN" + to_string(N_spins) + ".dat";
    
    string file_si_syntetic = "../tests/synthetic/si/si_synteticN" + to_string(N_spins) + ".dat";
    string file_sisj_syntetic = "../tests/synthetic/sisj/sisj_synteticN" + to_string(N_spins) + ".dat";
    string file_sisjsk_syntetic = "../tests/synthetic/sisjsk/sisjsk_synteticN" + to_string(N_spins) + ".dat";
    string file_Pij_syntetic = "../tests/synthetic/Pij/Pij_synteticN" + to_string(N_spins) + ".dat";
    string file_Cij_syntetic = "../tests/synthetic/Cij/Cij_synteticN" + to_string(N_spins) + ".dat";
    string file_Tijk_syntetic = "../tests/synthetic/Tijk/Tijk_synteticN" + to_string(N_spins) + ".dat";
    
    // Save synthetic data
    string file_states_syntetic = "../tests/data_synteticN" + to_string(N_spins) + ".dat";
    // Save synthetic mag_corr
    string file_mag_corr_syntetic = "../tests/mag_corr_synteticN" + to_string(N_spins) +  ".dat";
    //string energies_syntetic = "../tests/energiesN" + to_string(N_spins) +  ".dat";
    
    // Save Hamiltonian values (each row corresponds to a sigma state)
    saveValues(file_H_syntetic, "H", hamiltonianValues);
    // Save h_i values as a column
    saveValues(file_h_syntetic, "h_i", net.h);
    // Save upper triangular J_ij matrix
    saveValues(file_j_syntetic, "J_ij", net.J);;
    // Save sigma states (each row = one state)
    saveSigmaStates(file_states_syntetic, sigmaStates);
    
    // Save energies a long of time
    //SaveEnergy(energies_syntetic, energies);
    
    
    // SAVE PROPERTIES ===========
    // Save first moment (si/magnetization)
    saveS(file_si_syntetic,"si" ,si_synthetic);
    // Save Second moment (sisj)
    saveS(file_sisj_syntetic,"sisj", sisj_synthetic);
    // Save Third moment (sisjsk)
    saveS(file_sisjsk_syntetic,"sisjsk", sisjsk_synthetic);
    // Save Covariance (Cij)
    saveS(file_Cij_syntetic,"Cij", Cij_synthetic);
    // Save Correlation (Pij)
    saveS(file_Pij_syntetic,"Pij", Pij_synthetic);
    // Save Triplet (Tijk)
    saveS(file_Tijk_syntetic,"Tijk", Tijk_synthetic);

    computeMagCorr(sigmaStates, si_synthetic, sisj_synthetic, Cij_synthetic);
    saveMagCorr(file_mag_corr_syntetic, si_synthetic, sisj_synthetic, Cij_synthetic);

    // Saving information about synthetic data
    string info = "../tests/infoN" + to_string(N_spins) + ".dat";
    ofstream file_info (info.c_str());
    file_info << "N_spins = " << N_spins << "\tN_samples = " << M_states << "\tt_eq = " << to_string(t_eq) << "\trelx = " << to_string(relx) << "\trept = " << to_string(rept) << "\tseed = " << seed << endl;
    file_info.close();

    return 0;
}
