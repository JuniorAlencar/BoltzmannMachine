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

int main(int argc, char *argv[]) {

    // Gaussian Parameters ======================================================================================================= /
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
    
    // Opening syntetic data ======================================================================================================= /
    string file_states_syntetic = "../tests/data_syntetic.csv";
    ofstream file_data (file_states_syntetic.c_str());

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
    
    // Filenames to results of method ======================================================================================================= /
    string h_string = "../tests/" + method + "/hi/" + "file_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str;
    string J_string = "../tests/" + method + "/Jij/" + "file_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str;
    string H_string = "../tests/" + method + "/H/" + "file_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str;
    string si_string = "../tests/" + method + "/si/" + "file_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str;
    string sisj_string = "../tests/" + method + "/sisj/" + "file_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str;
    
    cout << "BMfinal ..." << endl;
	
	ifstream rede (file_rede_input.c_str());
	
	rede >> n;
	
	VecDoub av_s(n, 0.0), av_ss(n*(n-1)/2, 0.0), C(n*(n-1)/2, 0.0);
	
	//Rede r(tamanho, media, desvio, k, type = 0(random) 1(tree), H = 0(no field) 1(ConstField) -1(ConstField) 2(RandField))
	Rede r(n, mean, sigma, k, 0, 1);	
	
	for (int i = 0; i < r.nbonds; i++)
	{
		rede >> av_ss[i] >> C[i];
		
		if (i < r.n)
			rede >> av_s[i];	
		 
	}


    return 0;
}
