#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <fmt/core.h>
#include "./include/nr3.h"
#include "./include/network.h"
#include "./include/forwardmethod.h"
#include "./include/InverseMethod.h"
#include "./include/LUdcmp.h"

using namespace std;

int main(int argc, char *argv[]) {
	srand (time(NULL));
	
	// Gaussian Parameters ======================================================================================================= /
	int n;
	double mean = 1.0;
	double sigma = 0.25;
	double k = 10;
	int type;
	int H;
	
	int N_spins	= std::stoi(argv[1]);
	double min_erro_j	= std::stod(argv[2]);
	double min_erro_h	= std::stod(argv[3]);
	int multiply_teq 	= std::stoi(argv[4]);
	int multiply_relx 	= std::stoi(argv[5]);
    string method = argv[6];

    
    if (argc < 7) {
        std::cerr << "Uso: " << argv[0] << " <param1> <min_erro_j> <min_erro_h> <multi_teq> <multi_relx> <exact_solutions>" << std::endl;
        return 1;
    }

    try {
        double min_erro_j = std::stod(argv[2]);
        double min_erro_h = std::stod(argv[3]);
		int multiply_teq 	= std::stoi(argv[4]);
		int multiply_relx 	= std::stoi(argv[5]);

        cout << "min_erro_j: " << min_erro_j << endl;
        cout << "min_erro_h: " << min_erro_h << endl;
		cout << "multi_teq: " << multiply_teq << endl;
		cout << "multi_relx: " << multiply_relx << endl;
    } catch (const invalid_argument& e) {
        cerr << "Argumento inválido: " << e.what() << endl;
        return 1;
    } catch (const out_of_range& e) {
        cerr << "Valor fora do intervalo: " << e.what() << endl;
        return 1;
    }

    // Convert to string ======================================================================================================= /
	ostringstream os_j;
    os_j << scientific << setprecision(2) << min_erro_j;

	ostringstream os_h;
    os_h << scientific << setprecision(2) << min_erro_h;
	
	ostringstream os_teq;
	os_teq << multiply_teq;

	ostringstream os_relx;
	os_relx << multiply_relx;
	
	string min_erro_j_str = os_j.str();
	string min_erro_h_str = os_h.str();
	string multi_teq_str = os_teq.str();
	string multi_relx_str = os_relx.str();

	string errors_str = "../tests/" + method + "/erros_j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
    // Filenames to results of method ======================================================================================================= /
    // string h_string = "../tests/" + method + "/hi/" + "file_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str;
    // string J_string = "../tests/" + method + "/Jij/" + "file_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str;
    // string H_string = "../tests/" + method + "/H/" + "file_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str;
    // string si_string = "../tests/" + method + "/si/" + "file_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str;
    // string sisj_string = "../tests/" + method + "/sisj/" + "file_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str;
    
	
	cout << "Tests with synthetic data with N = " + to_string(N_spins) + " spins" << endl;
	
	string file_network_input = "../tests/mag_corr_synteticN" + to_string(N_spins) + ".dat";
	
	ifstream network (file_network_input.c_str());

	if (!network.is_open()) {
        std::cerr << "Erro to open file." << std::endl;
        return 1;
    }
	
	int N_pairs = (N_spins * (N_spins - 1)) / 2;
	int N_triplets = (N_spins * (N_spins - 1) * (N_spins - 2)) / 6;
	
	network >> n;
	
	// Vectors with <s>, <ss> and <C> from synthetic data ====================
	VecDoub av_s(N_spins, 0.0), av_ss(N_pairs, 0.0), C(N_pairs, 0.0);

	// Transfer values experimental to VecDoub
 	double ss, c, s;
    for (int i = 0; i < N_pairs; ++i) {
        if (!(network >> ss >> c >> s)) {
            std::cerr << "Error reading line " << i + 2 << std::endl;
            return 1;
        }

        av_ss[i] = ss;
        C[i] = c;

        if (i < N_spins) {
            av_s[i] = s;
        }
    }

    network.close();
	
	// Initial network to Boltmann Machine
	Rede bm(N_spins, 0, 0, 0, 0, 0);
	
	// <s> and <ss> to update with Ising
	VecDoub_IO bm_av_s(N_spins, 0.0), bm_av_ss(N_pairs, 0.0);

	// MC variables
	int t_eq = n*multiply_teq; // 150
	int relx = n*multiply_relx; // 2
	int rept = 40;
	int t_step = n*6000*relx/rept;
	
	double erroJ = 1, erroh = 1;
	double dJ, dh;
	int cort = 1000;
	int inter = 1;//= inter_ini;
	int inter_max = 300000;
	
	double eta_J = 0.05;//atof(argv[2]);
	double eta_h = 0.03;
	
	ofstream erros (errors_str.c_str());
	
	erros << "inter" << " " <<  "erroJ" << " " << "erroh" << endl; 
	
	
	while ((erroJ > min_erro_j || erroh > min_erro_h) && inter <= inter_max)    //(inter <= inter_max)
	{	
	srand(time(NULL)*time(NULL));
	erroJ = erroh = 0;

	eta_J = pow(inter, -0.4);
	eta_h = 2*pow(inter, -0.4);
	if(method == "metropolis")
		metropolis_bm(bm, bm_av_s, bm_av_ss, t_eq, t_step, relx, rept, 1);
	if(method == "exact" && N_spins < 25)
		exact_solution_bm (bm, bm_av_s, bm_av_ss, 1);

	for (int i = 0; i < bm.nbonds; i++)
	{
		if (i < bm.n)
		{
			dh = eta_h*(bm_av_s[i] - av_s[i]);
			erroh += pow(bm_av_s[i] - av_s[i], 2);
			bm.h[i] -= dh;
		}
		
		dJ = eta_J*(bm_av_ss[i] - av_ss[i]);
		erroJ += pow(bm_av_ss[i] - av_ss[i], 2);
		bm.J[i] -= dJ;

	}
	
	erroJ = sqrt(erroJ/bm.nbonds);
	erroh = sqrt(erroh/bm.n);

	//Salvando Erros
	erros << inter << " " << setprecision(13) << erroJ << " " << setprecision(13) << erroh << endl; 
	
	if (inter%cort == 0 || (erroJ < min_erro_j && erroh < min_erro_h))
	{
		std::cout << N_spins << " " << inter << " "
			<< "err_J" << " " << left << setw(13) << scientific << setprecision(6) << erroJ << " "
			<< "err_h" << " " << left << setw(13) << scientific << setprecision(6) << erroh << '\n';
	}			

	inter++;

	}

	erros.close();
    
	return 0;
}
