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
	
    string method = argv[1];
    int N_spins = stoi(argv[2]);
    int M_states = stoi(argv[3]);
    // double min_erro_j	= std::stod(argv[4]);
	// double min_erro_h	= std::stod(argv[5]);
	// int multiply_teq 	= std::stoi(argv[6]);
	// int multiply_relx 	= std::stoi(argv[7]);
    int seed = stoi(argv[4]);
    
    // Lista de métodos válidos
    vector<string> valid_methods = {
        "exact", "metropolis", "parallel_tempering", "swendsen_wang", "wang_landau"
    };

    // if (argc < 9) {
    //     cerr << "Uso: " << argv[0] << " <method> <N_spins> <M_states>" << endl;
    //     cerr << "Métodos disponíveis: exact, metropolis, parallel_tempering, swendsen_wang, wang_landau" << endl;
    //     return 1;
    // }

    // Verifica se o método fornecido é válido
    if (find(valid_methods.begin(), valid_methods.end(), method) == valid_methods.end()) {
        cerr << "Erro: método '" << method << "' inválido." << endl;
        cerr << "Métodos válidos: exact, metropolis, parallel_tempering, swendsen_wang, wang_landau" << endl;
        return 1;
    }
    
    cout << "Selected method: " << method << endl;
    cout << "Number of spins: " << N_spins << ", Number of states: " << M_states << endl;
    
    // Number of combinations in pairs s_i * s_j (s_i * s_j = s_j * s_i) -> symmetric matrix
    int N_pairs = (N_spins * (N_spins - 1)) / 2;
    
    // Initial network, with J and h obtaining from gaussian and uniform distributions, respectively, with mean=0.0
    Rede net_ini(N_spins, mean, sigma, k, 0, 1, seed);
    
    // av_s, av_ss is the initial avarage s_i and s_i*s_j
    vector<double> av_s = net_ini.h, av_ss = net_ini.J;
    
    // Network to update from MC method, using BoltzmannMachine(BM)
    Rede net_upd(N_spins, 0, 0, 0, 0, 0, seed);

    // bm_s, bm_ss is the update avarage s_i and s_i*s_j, using MC methods
    vector<double> bm_av_s(N_spins, 0.0), bm_av_ss(N_pairs, 0.0);
    
    // Vector to alocate the energies throughout the implementation
    vector<double> energies;
    vector<vector<int>> states;
    
    // MONTE CARLO UPDATE ----------------------
    int multiply_teq = 150;
    int multiply_relx = 2;
    // Monte Carlo variables
	int t_eq = N_spins*multiply_teq; // 150
	int relx = N_spins*multiply_relx; // 2
	int rept = 40;
	int t_step = n*6000*relx/rept;
    
	double min_error_j = 1.0e-5;
	double min_error_h = 1.0e-4;
    
    // Initial errors in J and h
    double erroJ = 1.0, erroh = 1.0;
    double dJ, dh;
	int cort = 200;
	// inter_max - inter = Number of interations
    int inter = 1;
	int inter_max = 300000;
	
    // Process update rate
	double eta_J = 0.05;
	double eta_h = 0.03;
    string file_name_errors = "../tests/errors.dat";
    string file_name_energies = "../tests/energies.dat";
    
    ofstream erros (file_name_errors.c_str());
    ofstream ener (file_name_energies.c_str());
    
    erros << "inter" << " " <<  "erroJ" << " " << "erroh" << endl; 
    ener << "inter" << " " << "energy" << endl;
    
    while ((erroJ > min_error_j || erroh > min_error_h) && inter <= inter_max) {
        srand(time(NULL) * time(NULL));

        erroJ = erroh = 0;

        eta_J = pow(inter, -0.4);
        eta_h = 2 * pow(inter, -0.4);

        GenerateStates(net_upd, bm_av_s, bm_av_ss, t_eq,  relx, rept, M_states, 1.0, energies, states, seed);

        for (int i = 0; i < net_upd.nbonds; i++) {
            if (i < net_upd.n) {
                dh = eta_h * bm_av_s[i];
                erroh += pow(bm_av_s[i], 2);
                net_upd.h[i] -= dh;
            }
            dJ = eta_J*bm_av_ss[i];
			erroJ += pow(bm_av_ss[i], 2);
			net_upd.J[i] -= dJ;
        }
        erroJ = sqrt(erroJ / net_upd.nbonds);
        erroh = sqrt(erroh / net_upd.n);

        erros << inter << " " << setprecision(13) << erroJ << " " << erroh << endl;
        
        // Média da energia dessa iteração
        double ava_energy = 0.0;
        for (size_t i = 0; i < energies.size(); ++i)
            ava_energy += energies[i];
        ava_energy /= energies.size();
        
        // Salvando energia no mesmo estilo
        ener << inter << " " << setprecision(13) << ava_energy << endl;        
        
        // if (inter % cort == 0 || (erroJ < min_error_j && erroh < min_error_h)) {
        //     cout << "tests " << inter << " "
        //          << "err_J " << left << setw(13) << scientific << setprecision(6) << erroJ << " "
        //          << "err_h " << left << setw(13) << scientific << setprecision(6) << erroh << endl;
        // }

        inter++;
    }

    erros.close();
    ener.close();
    //for (int i = 0; i<N_spins;i++)
        //cout << h[i] << endl;
    // Rede r(N_spins, mean, sigma, k, 0, 1);
    // for (int i =0; i < r.nbonds; i++){
    //     J[i] >> r.J[i];
    //     if (i < N_spins)
    //         h[i] >> r.h[i];
    // }
    // for (int i =0; i < r.nbonds; i++){
    //     cout << r.J[i] << endl;
    //     if (i < N_spins)
    //         cout << r.h[i] << endl;
    // }
    // // Generate M sigma states without correlation
    // vector<double> hamiltonianValues = computeHamiltonian(sigmaStates, h, J);

    // // Calculate si and sisj and C_i from synthetic data
    // vector<double> si_synthetic = computeSi(sigmaStates);
    // vector<double> sisj_synthetic = computeSiSj(sigmaStates);
    // vector<double> C_synthetic = computeC(si_synthetic, sisj_synthetic);
    
    // string file_h_syntetic = "../tests/" + method +  "/hi/h_syntetic.csv";
    // string file_si_syntetic = "../tests/" + method +  "/si/si_syntetic.csv";
    // string file_sisj_syntetic = "../tests/" + method +  "/sisj/sisj_syntetic.csv";
    // string file_j_syntetic = "../tests/" + method +  "/Jij/J_syntetic.csv";
    // string file_H_syntetic = "../tests/" + method +  "/H/H_syntetic.csv";
    
    // // Save synthetic data
    // string file_states_syntetic = "../tests/data_synteticN" + to_string(N_spins) + ".dat";
    // // Save synthetic mag_corr
    // string file_mag_corr_syntetic = "../tests/mag_corr_synteticN" + to_string(N_spins) +  ".dat";
    
    // // Save h_i values as a column
    // savehH(file_h_syntetic, h);

    // // Save upper triangular J_ij matrix
    // saveJ(file_j_syntetic, J);

    // // Save sigma states (each row = one state)
    // saveSigmaStates(file_states_syntetic, sigmaStates);

    // // Save Hamiltonian values (each row corresponds to a sigma state)
    // savehH(file_H_syntetic, hamiltonianValues);

    // saveSi(file_si_syntetic, si_synthetic);
    // saveSiSj(file_sisj_syntetic, sisj_synthetic);
    // //saveMagCorr(file_mag_corr_syntetic, si_synthetic, sisj_synthetic, C_synthetic);
    
    // computeMagCorr(sigmaStates, si_synthetic, sisj_synthetic, C_synthetic);
    // saveMagCorr(file_mag_corr_syntetic, si_synthetic, sisj_synthetic, C_synthetic);

    // // Saving information about synthetic data
    // string info = "../tests/" + method + "/info.dat";
    // ofstream file_info (info.c_str());
    // file_info << "N_spins = " << N_spins << "\tN_samples = " << M_states << endl;
    // file_info.close();

    return 0;
}
