#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <cstdio>  // remove file
#include <fstream>
#include <fmt/core.h>
#include "./include/nr3.h"
#include "./include/network.h"
#include "./include/forwardmethod.h"
#include "./include/InverseMethod.h"
#include "./include/LUdcmp.h"
#include "./include/create_folders.h"
#include "./include/exp_means.h" // Calculate experimental means
#include "./include/write_json.h" // Write results in json

using namespace std;

int main(int argc, char *argv[]){
	
	srand (time(NULL));
	
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
	bool use_exact = (std::string(argv[6]) == "true");
	
    if (argc < 7) {
        cerr << "Uso: " << argv[0] << " <param1> <min_erro_j> <min_erro_h> <multi_teq> <multi_relx> <exact_solutions>" << endl;
        return 1;
    }

   	// Load experimental means from data
	struct exp_means My_exp = experimental_means(text_name, multiply_teq, multiply_relx, min_erro_h, min_erro_j, use_exact);
	
	// Converter para string
	ostringstream os_j;
    os_j << scientific << setprecision(2) << min_erro_j;

	ostringstream os_h;
    os_h << scientific << setprecision(2) << min_erro_h;
	
	ostringstream os_teq;
	os_teq << My_exp.t_eq;

	ostringstream os_relx;
	os_relx << My_exp.relx;

	string min_erro_j_str = os_j.str();
	string min_erro_h_str = os_h.str();
	string teq_str = os_teq.str();
	string relx_str = os_relx.str();
	
	// If true -> exact_solutions, else -> Metropolis--------------------
	string file_rede_input;
	string file_network_name;
	string file_rede_output;
	string file_name_erros;
	
	if(use_exact == true){
		file_rede_input = "../Data/Mag_Corr/" + text_name + "/exact/teq_" + teq_str + "/relx_" + relx_str + "/exp_j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + ".dat";
		file_network_name = "../Results/" + text_name + "/teq_" + teq_str + "/relx_" + relx_str + "/Network/j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + ".dat";
		file_rede_output = "../Results/" + text_name + "/teq_" + teq_str + "/relx_" + relx_str + "/Network/j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + ".dat";
		file_name_erros = "../Results/" + text_name + "teq_" + teq_str + "/relx_" + relx_str + "/Errors/j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + ".dat";

	}
	else{
		file_rede_input = "../Data/Mag_Corr/" + text_name + "/metropolis/teq_" + teq_str + "/relx_" + relx_str + "/exp_j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + ".dat";
		file_network_name = "../Results_Metropolis/" + text_name + "/teq_" + teq_str + "/relx_" + relx_str + "/Network/j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + ".dat";
		file_rede_output = "../Results_Metropolis/" + text_name + "/teq_" + teq_str + "/relx_" + relx_str + "/Network/j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + ".dat";
		file_name_erros = "../Results_Metropolis/" + text_name + "teq_" + teq_str + "/relx_" + relx_str + "/Errors/j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + ".dat";
	}
	// nomes dos arquivos

	//-----------------------------------------------------------------------------
	//Ler o arquivo com as correlações e magnetizações de um certo arquivo

	//create_folders();
	
	cout << "BMfinal ..." << endl;
	
	ifstream rede (file_rede_input.c_str());
	
	rede >> n;
	// First moment, second moment
	VecDoub_IO bm_s(n, 0.0), bm_ss(n*(n-1)/2, 0.0), C(n*(n-1)/2, 0.0);
	
	//Rede r(tamanho, media, desvio, k, type = 0(random) 1(tree), H = 0(no field) 1(ConstField) -1(ConstField) 2(RandField))
	Rede r(n, mean, sigma, k, 0, 1);	
	
	for (int i = 0; i < r.nbonds; i++)
	{
		rede >> bm_ss[i] >> C[i];
		
		if (i < r.n)
			rede >> bm_s[i];	
		 
	}
	
	//rede.close();
	
//-----------------------------------------------------------------------------
//Abrir arquivo da rede

	ifstream network_in (file_network_name.c_str());
	
	network_in >> r.n;
	
	//rede para ser atualizada
	Rede bm(n, 0, 0, 0, 0, 0);
	
	for(int i = 0; i < r.nbonds; i++)
	{
		network_in >> bm.J[i];

		if (i < n)
			network_in >> bm.h[i];
	
	}
	
	network_in.close();

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

	//rede para ser atualizada
	//Rede bm(n, mean, sigma, k, 0, 1);

	VecDoub_IO am_av_s(n, 0.0), am_av_ss(n*(n-1)/2, 0.0);

	//variaveis para MC
	int t_eq = n*multiply_teq; // 300
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

	//Arquivo para salvar os erros ao longo do tempo
	
	ofstream erros (file_name_erros.c_str());
	
	erros << "inter" << " " <<  "erroJ" << " " << "erroh" << endl; 
	
	while ((erroJ > min_erro_j || erroh > min_erro_h) && inter <= inter_max)    //(inter <= inter_max)
	{	
		srand(time(NULL)*time(NULL));

		erroJ = erroh = 0;

		eta_J = pow(inter, -0.4);
		eta_h = 2*pow(inter, -0.4);

		
		if (use_exact == false)
		{
			metropolis_bm (bm, am_av_s, am_av_ss, t_eq, t_step, relx, rept, 1);
		}
		else
		{
			exact_solution_bm (bm, am_av_s, am_av_ss, 1);
		}		

		for (int i = 0; i < bm.nbonds; i++)
		{
			if (i < bm.n)
			{
				dh = eta_h*(am_av_s[i] - bm_s[i]);
				erroh += pow(am_av_s[i] - bm_s[i], 2);
				bm.h[i] -= dh;
			}
			
			dJ = eta_J*(am_av_ss[i] - bm_ss[i]);
			erroJ += pow(am_av_ss[i] - bm_ss[i], 2);
			bm.J[i] -= dJ;

		}
		
		erroJ = sqrt(erroJ/bm.nbonds);
		erroh = sqrt(erroh/bm.n);

		//Salvando Erros
		erros << inter << " " << setprecision(13) << erroJ << " " << setprecision(13) << erroh << endl; 
		
		if (inter%cort == 0 || (erroJ < min_erro_j && erroh < min_erro_h))
		{
			std::cout << text_name << " " << inter << " "
              << "err_J" << " " << left << setw(13) << scientific << setprecision(6) << erroJ << " "
              << "err_h" << " " << left << setw(13) << scientific << setprecision(6) << erroh << '\n';
		}			

		inter++;
	
	}
	
	//Fechar arquivos dos erros salvos 
	erros.close();
	
//-----------------------------------------------------------------------------
	//Arquivo para salvar a rede obtida	
	ofstream network (file_rede_output.c_str());

	network << bm.n << endl;
	
	for (int i = 0; i < bm.nbonds; i++)
	{
		network << left << setw(15) << bm.J[i];
		
		if (i < bm.n)
			network << left << setw(15) << bm.h[i];
			
		network << endl;
	}
	
	network.close();
	
//-----------------------------------------------------------------------------
//Arquivo para salvar a correlação e magnetizações geradas pela rede encontrada

	// // Remove file Mag_Corr
    // if (remove(file_rede_input.c_str()) == 0) {
    //     std::cout << "Mag_Corr_Exp delete!" << std::endl;
    // } else {
    //     std::perror("Erro to delete Mag_Corr_Exp");
    // }
	
	VecDoub am_C(n*(n-1)/2, 0.0);

	//correlação pearson
	vector<double> pearson_ising(n*(n-1)/2);
	
	int ind = 0;
	for (int p = 0; p < n-1; p++)
	{
		for (int pp = p+1; pp < n; pp++)
		{
			am_C[ind] = am_av_ss[ind] - am_av_s[p]*am_av_s[pp];

			pearson_ising[ind] = am_C[ind]/sqrt((1 - pow(am_av_s[p], 2))*(1 - pow(am_av_s[pp], 2)));
			
			ind++;
		}
	}
		
//-----------------------------------------------------------------------------
//Salvar os dados separadamente

	//-----------------------------------------------------
	//Calcular os tripletos

	int n_triplet = n*(n-1)*(n-2)/6;
	vector<double> am_av_sss(n_triplet), ising_Triplet(n_triplet);
	vector<vector<double>> ising_M_av_ss(n, vector<double>(n));

	// Resolvendo dependendo do numero de observaveis
	if (use_exact == false)
	{
		metropolis_triplet (bm, am_av_s, am_av_ss, am_av_sss, t_eq, t_step, relx, rept, 1);
	}
	else
	{
		exact_solution_triplet (bm, am_av_s, am_av_ss, am_av_sss, 1);
	}


	//Interação par-a-par 
	ind = 0;
	for (int p = 0; p < n-1; p++)
	{
		for(int pp = p+1; pp < n; pp++)
		{		

			ising_M_av_ss[p][pp] = am_av_ss[ind];
			ising_M_av_ss[pp][p] = ising_M_av_ss[p][pp];
		
			ind++;
		}
	}


	//Terceiro momento centrado (tripleto)
	ind = 0;
	for (int i = 0; i < n-2; i++)
	{
		for (int j = i+1; j < n-1; j++)
		{
			for (int k = j+1; k < n; k++)
			{
				ising_Triplet[ind] = am_av_sss[ind] - am_av_s[i]*ising_M_av_ss[j][k] - am_av_s[j]*ising_M_av_ss[i][k] - am_av_s[k]*ising_M_av_ss[i][j]
									 + 2*am_av_s[i]*am_av_s[j]*am_av_s[k];
								
				ind++;
			}
		}
	}

	
	// Convert VecDoub to vector<double> (first and second moments)
	vector<double> am_s(bm.n, 0), am_ss(bm.nbonds,0), amC(n*(n-1)/2);
	for (int i = 0; i < bm.n; ++i)
		am_s[i] = am_av_s[i];
	for (int i = 0; i < bm.nbonds; ++i)
		am_ss[i] = am_av_ss[i];
	for (int i = 0; i < n*(n-1)/2; ++i)
		amC[i] = am_C[i];
	
	// Write json file with results
	write_json_properties(text_name, relx, t_eq, My_exp.nspins, use_exact, 
						  min_erro_h, min_erro_j, bm ,My_exp, am_av_s, am_av_ss, am_av_sss, am_C, pearson_ising,
						  ising_Triplet);
	//My_exp
	
	return 0;
}
