#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include "./include/nr3.h"
#include "./include/network.h"
#include "./include/forwardmethod.h"
#include "./include/InverseMethod.h"
#include "./include/LUdcmp.h"

using namespace std;

int main(int argc, char *argv[]){

	srand (time(NULL));

	int n;
	double mean = 1;
	double sigma = 0.25;
	double k = 10;
	int type;
	int H;
	
	//string text_input  = argv[1];
	//string rede_output = argv[2];
	string text_name     = argv[1];
	
//-----------------------------------------------------------------------------
//Ler o arquivo com as correlações e magnetizações de um certo arquivo

	string file_rede_input = "../Data/Mag_Corr/mag_corr_exp_" + text_name + ".dat";

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
	
	rede.close();
	
//-----------------------------------------------------------------------------
//Verifica os vetores de magnetizações e correlações

//	for (int i = 0; i < r.nbonds; i++)
//	{
//		cout << left << setw(15) << av_ss[i] << left << setw(15) << C[i];
//		
//		if (i < r.n)
//			cout << left << setw(15) << av_s[i];
//			
//		cout << endl;
//	}

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//Abrir arquivo da rede

	string file_network_name = "../Results/Network/network_" + text_name + ".dat";
	
	ifstream network_in (file_network_name.c_str());
	
	network_in >> n;

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

	VecDoub bm_av_s(n, 0.0), bm_av_ss(n*(n-1)/2, 0.0);

	//variaveis para MC
	int t_eq = n*150;
	int relx = 2*n;
	int rept = 40;
	int t_step = n*6000*relx/rept;
	
	double erroJ = 1, erroh = 1;
	double dJ, dh;
	//int inter_ini = atoi(argv[3]);
	//int inter_max = atoi(argv[4]);
	int cort = 1000;
	int inter = 1;//= inter_ini;
	int inter_max = 100000;
	
	double eta_J = 0.05;//atof(argv[2]);
	double eta_h = 0.03;

	// N30 (j,h) -> (4e-5, 2e-4)
	// N30 new (j,h) -> (9e-6, 8e-5)
	
	// N20 (j,h) -> (6e-7, 5e-6)
	double min_erro_j = 2.0e-7;
	double min_erro_h = 2.0e-6;

	bool use_exact = false;

	//Arquivo para salvar os erros ao longo do tempo
	string file_name_erros = "../Results/Erro/erro_" + text_name + ".dat";

	ofstream erros (file_name_erros.c_str());

	//erros.seekg(0, std::ios_base::end);
	
	while ((erroJ > min_erro_j || erroh > min_erro_h) && inter <= inter_max)    //(inter <= inter_max)
	{	
		srand(time(NULL)*time(NULL));

		erroJ = erroh = 0;

		eta_J = pow(inter, -0.4);
		eta_h = 2*pow(inter, -0.4);

		
		if (n > 25 && use_exact == false)
		{
			metropolis_bm (bm, bm_av_s, bm_av_ss, t_eq, t_step, relx, rept, 1);
		}
		else
		{
			exact_solution_bm (bm, bm_av_s, bm_av_ss, 1);
		}		

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
		
		erroJ = sqrt(erroJ)/bm.nbonds;
		erroh = sqrt(erroh)/bm.n;

		//Salvando Erros
		erros << inter << " " << setprecision(13) << erroJ << " " << setprecision(13) << erroh << endl; 
		
		if (inter%cort == 0 || (erroJ < min_erro_j && erroh < min_erro_h))
		{
			cout  << text_name << " " << inter << "  " << left << setw(13) << erroJ << left << setw(13) << erroh << endl;
			//save_data (r, bm, av_s, av_ss, bm_av_s, bm_av_ss, pasta, inter);

			//salvando dados
			//SalvarDados(bm, bm_av_s, bm_av_ss, text_name);
		}			

		inter++;
	
	}
	
	//Fechar arquivos dos erros salvos 
	erros.close();
	
//-----------------------------------------------------------------------------
//Arquivo para salvar a rede obtida	

	string file_rede_output = "../Results/Network/network_" + text_name + ".dat";

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

	string file_mag_corr_output = "../Results/Mag_Corr-ising/mag_corr_ising_" + text_name + ".dat";

	ofstream mag_corr (file_mag_corr_output.c_str());

	VecDoub bm_C(n*(n-1)/2, 0.0);

	//correlação pearson
	vector<double> pearson_ising(n*(n-1)/2);
	
	int ind = 0;
	for (int p = 0; p < n-1; p++)
	{
		for (int pp = p+1; pp < n; pp++)
		{
			bm_C[ind] = bm_av_ss[ind] - bm_av_s[p]*bm_av_s[pp];

			pearson_ising[ind] = bm_C[ind]/sqrt((1 - pow(bm_av_s[p], 2))*(1 - pow(bm_av_s[pp], 2)));
			
			ind++;
		}
	}

	mag_corr << bm.n << endl;
	
	for (int i = 0; i < bm.nbonds; i++)
	{
		mag_corr << left << setw(15) << bm_av_ss[i] << left << setw(15) << bm_C[i];
		
		if (i < bm.n)
			mag_corr << left << setw(15) << bm_av_s[i];
			
		mag_corr << endl;
	}
	
	mag_corr.close();

	//Salva arquivo com Jij e Correlação-------------------
    //Nome do arquivo alvo
    string file_name_CorrJij = "../Results/CorrJij/CorrJij_" + text_name + ".dat";

    //Abrindo arquivo output
    ofstream CorrJij (file_name_CorrJij.c_str());

	//Salvar arquivo com Jij e Pij
	//Nome do arquivo alvo
	string file_name_PJij = "../Results/PJij/PJij_" + text_name + ".dat";

	//Abrindo arquivo output
	ofstream PJij (file_name_PJij.c_str());

    //Passando os valores para o arquivo
    for (int i = 0; i < bm.nbonds; i++)
    {
        CorrJij << left << setw(20) << setprecision(13) << bm.J[i] << 
                left << setw(18) << setprecision(13) << bm_C[i] << endl;

		PJij << bm.J[i] << " " << pearson_ising[i] << endl;
    }
	
//-----------------------------------------------------------------------------
//Salvar os dados separadamente

    //Arquivo para hi------------------
    string file_name_hi = "../Results/SeparateData/hi/hi_" + text_name + ".dat";

    //abrir arquivo
    ofstream file_hi (file_name_hi.c_str());

    //Arquivo para mi -----------------
    string file_name_mi = "../Results/SeparateData/mi-ising/mi_ising_" + text_name + ".dat";

    //abrir arquivo
    ofstream file_mi (file_name_mi.c_str());

    for (int i = 0; i < n; i++)
    {
        file_hi << bm.h[i] << endl;
        file_mi << bm_av_s[i] << endl;
    }

    //fechando arquivos
    file_hi.close();
    file_mi.close();

    //-----------------------------------------------------

    //Arquivo para Jij------------------
    string file_name_Jij = "../Results/SeparateData/Jij/Jij_" + text_name + ".dat";

    //abrir arquivo
    ofstream file_Jij (file_name_Jij.c_str());

    //Arquivo para Cij -----------------
    string file_name_Cij = "../Results/SeparateData/Cij-ising/Cij_ising_" + text_name + ".dat";

    //abrir arquivo
    ofstream file_Cij (file_name_Cij.c_str());

	//Arquivo para Pij
	string file_name_Pij = "../Results/SeparateData/Pij-ising/Pij_ising_" + text_name + ".dat";

	//abrir arquivo
	ofstream file_Pij (file_name_Pij.c_str());

	//Arquivo para sisj
	string file_name_sisj = "../Results/SeparateData/sisj-ising/sisj_ising_" + text_name + ".dat";

	//abrir arquivo
	ofstream file_sisj (file_name_sisj.c_str());

    //salvar arquivos
    for (int i = 0; i < (n*(n-1)/2); i++)
    {
        file_Jij  << bm.J[i] << endl;
        file_Cij  << bm_C[i] << endl;
		file_Pij  << pearson_ising[i] << endl;
		file_sisj << bm_av_ss[i] << endl;
    }

    //fechando arquivos
    file_Jij.close();
    file_Cij.close();
	file_Pij.close();
    file_sisj.close();

	//-----------------------------------------------------
	//Calcular os tripletos

	int n_triplet = n*(n-1)*(n-2)/6;
	vector<double> bm_av_sss(n_triplet), ising_Triplet(n_triplet);
	vector<vector<double>> ising_M_av_ss(n, vector<double>(n));

	//Arquivo para Tijk
	string file_name_Tijk = "../Results/SeparateData/Tijk-ising/Tijk_ising_" + text_name + ".dat";

	//abrir arquivo
	ofstream file_Tijk (file_name_Tijk.c_str());

	//Arquivo para sisjsk
	string file_name_sisjsk = "../Results/SeparateData/sisjsk-ising/sisjsk_ising_" + text_name + ".dat";

	//abrir arquivo
	ofstream file_sisjsk (file_name_sisjsk.c_str());


	// Resolvendo dependendo do numero de observaveis
	if (n > 25)
	{
		metropolis_triplet (bm, bm_av_s, bm_av_ss, bm_av_sss, t_eq, t_step, relx, rept, 1);
	}
	else
	{
		exact_solution_triplet (bm, bm_av_s, bm_av_ss, bm_av_sss, 1);
	}


	//Interação par-a-par 
	ind = 0;
	for (int p = 0; p < n-1; p++)
	{
		for(int pp = p+1; pp < n; pp++)
		{		

			ising_M_av_ss[p][pp] = bm_av_ss[ind];
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
				ising_Triplet[ind] = bm_av_sss[ind] - bm_av_s[i]*ising_M_av_ss[j][k] - bm_av_s[j]*ising_M_av_ss[i][k] - bm_av_s[k]*ising_M_av_ss[i][j]
									 + 2*bm_av_s[i]*bm_av_s[j]*bm_av_s[k];
				
				file_Tijk   << ising_Triplet[ind] << endl;
				file_sisjsk << bm_av_sss[ind] << endl;
				
				ind++;
			}
		}
	}

	
	file_Tijk.close();
    file_sisjsk.close();

	return 0;
}
