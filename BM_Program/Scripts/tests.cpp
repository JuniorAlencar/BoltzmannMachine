#include "include/synthetic_H.h"

int main(int argc, char *argv[]) {

    // Gaussian Parameters
	int n;
	double mean = 1;
	double sigma = 0.25;
	double k = 10;
	int type;
	int H;
	
	string text_name	= argv[1];
	double min_erro_j	= std::stod(argv[2]);
	double min_erro_h	= std::stod(argv[3]);
	int multiply_teq 	= std::stoi(argv[4]);
	int multiply_relx 	= std::stoi(argv[5]);
    string method = argv[6];
    
    int N_spins = 30; // Number of spins
    int M_states = 50; // Number of sigma states
    
    // Generate h_i values
    vector<double> h = generateHVector(N_spins, -1.0, 1.0);

    // Generate J_ij values
    vector<vector<double>> J = generateJMatrix(N_spins, mean, sigma);

    // Generate M sigma states
    vector<vector<double>> sigmaStates;
    vector<double> hamiltonianValues;

    for (int m = 0; m < M_states; ++m) {
        vector<double> sigma = generateBinaryVector(N_spins);
        sigmaStates.push_back(sigma);

        // Compute Hamiltonian for this state
        double H = computeHamiltonian(sigma, h, J);
        hamiltonianValues.push_back(H);
    }

    // Calculate si and sisj from synthetic data
    vector<double> si_synthetic = computeSi(sigmaStates);
    vector<double> sisj_synthetic = computeSiSj(sigmaStates);
    
    string file_h_syntetic = "../tests/" + method +  "/hi/h_syntetic.csv";
    string file_si_syntetic = "../tests/" + method +  "/si/si_syntetic.csv";
    string file_sisj_syntetic = "../tests/" + method +  "/sisj/sisj_syntetic.csv";
    string file_j_syntetic = "../tests/" + method +  "/Jij/J_syntetic.csv";
    string file_H_syntetic = "../tests/" + method +  "/H/H_syntetic.csv";
    string file_states_syntetic = "../tests/data_syntetic.csv";
    
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


    return 0;
}
