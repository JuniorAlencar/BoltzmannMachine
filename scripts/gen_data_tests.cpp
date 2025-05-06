#include "./include/synthetic_data.h"

#include <algorithm>

int main(int argc, char *argv[]) {

    // Gaussian Parameters
	int n;
	double mean = 0.0;
	double sigma = 0.7; // N = 60
    //double sigma = 0.1; //N = 20, N = 40
	double k = 10;
	int type;
	
    int N_spins = stoi(argv[1]);
    int M_states = stoi(argv[2]);
    string method = argv[3];
    int seed = stoi(argv[4]);
    int multi_relx = stoi(argv[5]);
    int multi_teq = stoi(argv[6]);
    string tests = argv[7]; // On or off
    
    if (argc < 7) {
        cerr << "Uso: " << argv[0] << "<N_spins> <M_states> <seed> <tests>" << endl;
        
        return 1;
    }
    cout << "Number of spins: " << N_spins << ", Number of states: " << M_states << ", seed: " << seed << endl;
    
    // Generate to code
    std::mt19937 gen(seed);
    
    // Number of combinations in pairs s_i * s_j (s_i * s_j = s_j * s_i) -> symmetric matrix
    int N_pairs = (N_spins * (N_spins - 1)) / 2;
    
    // Initial network, with J and h obtaining from gaussian and uniform distributions, respectively, with mean=0.0
    Rede net(N_spins, mean, sigma, k, type, 1, gen);
    double min = -1.0, max = 1.0;
    
    vector<double> J_ij(N_pairs, 0.0);
    vector<double> h_i(N_spins, 0.0);
    
    for (int i = 0; i<N_pairs; i++){
        J_ij[i] = net.J[i];
        if(i < N_spins)
            h_i[i]= net.h[i];
    }
    
    vector<double> av_s(N_spins, 0.0), av_ss(N_pairs, 0.0);
    
    // Vector to alocate the energies throughout the implementation
    vector<double> energies;
    vector<vector<int>> sigmaStates;
    
    // MONTE CARLO UPDATE ----------------------
    
    // Monte Carlo variables
	int t_eq = multi_teq*N_spins;
	int relx = multi_relx*N_spins;
	int rept = 50;
    // N = 20 -> relx = 5, t_eq = 50
    // N = 40 -> relx = 10, t_eq = 50
    // N = 60 -> relx = 40, t_eq = 100
    // N = 80 -> relx = 50, t_eq = 150
    // Number of steps
	int t_step = (M_states * relx) / rept;
    
    GenerateStates_Wolff(net, av_s, av_ss, t_eq,  relx, rept, M_states, 1.0, energies, sigmaStates, gen);
    
    // Compute Hamiltonian
    vector<double> hamiltonianValues = computeHamiltonian(sigmaStates, h_i, J_ij);

    // Calculate si and sisj and C_i from synthetic data
    vector<double> si_synthetic = computeSi(sigmaStates);
    vector<double> sisj_synthetic = computeSiSj(sigmaStates);
    vector<double> Cij_synthetic = computeC(si_synthetic, sisj_synthetic);
    vector<double> Pij_synthetic = computePij(sigmaStates, Cij_synthetic);
    vector<double> sisjsk_synthetic = computeSiSjSk(sigmaStates);
    vector<double> Tijk_synthetic = computeTriplet(sigmaStates, sisjsk_synthetic);

    // Declate the folder results
    string file_h_syntetic;
    string file_j_syntetic;
    string file_si_syntetic;
    string file_sisj_syntetic;
    string file_sisjsk_syntetic;
    string file_Pij_syntetic;
    string file_Cij_syntetic;
    string file_Tijk_syntetic;

    string file_states_syntetic;
    string file_mag_corr_syntetic;

    string info;
    
    if(tests == "off"){
        file_h_syntetic = "../Results/" + method + "/SeparateData/hi/hi_real_data_synteticN" + to_string(N_spins) + ".dat";
        //file_H_syntetic = "../Results/" + method + "/H/H_synteticN/"
        file_j_syntetic = "../Results/" + method + "/SeparateData/Jij/Jij_real_data_synteticN" + to_string(N_spins) + ".dat";
        
        file_si_syntetic = "../Results/" + method + "/SeparateData/mi-exp/mi_exp_data_synteticN" + to_string(N_spins) + ".dat";
        file_sisj_syntetic = "../Results/" + method + "/SeparateData/sisj-exp/sisj_exp_data_synteticN" + to_string(N_spins) + ".dat";
        file_sisjsk_syntetic = "../Results/" + method + "/SeparateData/sisjsk-exp/sisjsk_exp_data_synteticN" + to_string(N_spins) + ".dat";
        file_Pij_syntetic = "../Results/" + method + "/SeparateData/Pij-exp/Pij_exp_data_synteticN" + to_string(N_spins) + ".dat";
        file_Cij_syntetic = "../Results/" + method + "/SeparateData/Cij-exp/Cij_exp_data_synteticN" + to_string(N_spins) + ".dat";
        file_Tijk_syntetic = "../Results/" + method + "/SeparateData/Tijk-exp/Tijk_exp_data_synteticN" + to_string(N_spins) + ".dat";
        
        // Save synthetic data
        file_states_syntetic = "../Data/TidyData/data_synteticN" + to_string(N_spins) + ".dat";
        // Save synthetic mag_corr
        file_mag_corr_syntetic = "../Data/Mag_Corr/mag_corr_exp_data_synteticN" + to_string(N_spins) +  ".dat";
        // SAVE PROPERTIES ===========
        // Save first moment (si/magnetization)
        saveS(file_si_syntetic,"si_synt" ,si_synthetic);
        // Save Second moment (sisj)
        saveS(file_sisj_syntetic,"sisj_synt", sisj_synthetic);
        // Save Third moment (sisjsk)
        saveS(file_sisjsk_syntetic,"sisjsk_synt", sisjsk_synthetic);
        // Save Covariance (Cij)
        saveS(file_Cij_syntetic,"Cij_synt", Cij_synthetic);
        // Save Correlation (Pij)
        saveS(file_Pij_syntetic,"Pij_synt", Pij_synthetic);
        // Save Triplet (Tijk)
        saveS(file_Tijk_syntetic,"Tijk_synt", Tijk_synthetic);

        computeMagCorr(sigmaStates, si_synthetic, sisj_synthetic, Cij_synthetic);
        saveMagCorr(file_mag_corr_syntetic, si_synthetic, sisj_synthetic, Cij_synthetic);
        
        info = "../Data/TidyData/info_synteticN" + to_string(N_spins) + ".dat";
    }
    
    else if (tests== "on")
    {
        file_h_syntetic = "../syntetic_tests/N_" + to_string(N_spins) + "/hi/hi_" + to_string(seed) + ".dat";
        file_j_syntetic = "../syntetic_tests/N_" + to_string(N_spins) + "/Jij/Jij_" + to_string(seed) + ".dat";
        file_states_syntetic = "../syntetic_tests/N_" + to_string(N_spins) + "/data_" + to_string(seed) + ".dat";

        // Saving information about synthetic data
        info = "../syntetic_tests/N_" + to_string(N_spins) + "/Info/info_" + to_string(seed) + ".dat";

    }

    // Save Hamiltonian values (each row corresponds to a sigma state)
    //saveValues(file_H_syntetic, "H", hamiltonianValues);
    // Save h_i values as a column
    saveValues(file_h_syntetic, "hi_synt", h_i);
    // Save upper triangular J_ij matrix
    saveValues(file_j_syntetic, "Jij_synt", J_ij);
    // Save sigma states (each row = one state)
    saveSigmaStates(file_states_syntetic, sigmaStates);
    
    // Save energies a long of time
    //SaveEnergy(energies_syntetic, energies);
    
    // Saving Info about the implementation
    ofstream file_info (info.c_str());
    file_info << "N_spins = " << N_spins << "\tN_samples = " << M_states << "\tt_eq = " << to_string(t_eq) << "\trelx = " << to_string(relx) << "\trept = " << to_string(rept) << "\tseed = " << seed << endl;
    file_info.close();




    return 0;
}