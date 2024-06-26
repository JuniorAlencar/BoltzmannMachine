#include "network.h"
#include "forwardmethod.h"
#include <vector>
using namespace std;

void SalvarDados (Rede bm, VecDoub_I av_s, VecDoub_I av_ss, string text_name)
{
	int n = bm.n;

	VecDoub bm_av_s = av_s;
	VecDoub bm_av_ss = av_ss;

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
	VecDoub_IO  bm_av_sss(n_triplet), ising_Triplet(n_triplet);
	vector<vector<double>> ising_M_av_ss(n, vector<double>(n));

	//Arquivo para Tijk
	string file_name_Tijk = "../Results/SeparateData/Tijk-ising/Tijk_ising_" + text_name + ".dat";

	//abrir arquivo
	ofstream file_Tijk (file_name_Tijk.c_str());

	//Arquivo para sisjsk
	string file_name_sisjsk = "../Results/SeparateData/sisjsk-ising/sisjsk_ising_" + text_name + ".dat";

	//abrir arquivo
	ofstream file_sisjsk (file_name_sisjsk.c_str());

	//variaveis para MC
	int t_eq = n*150;
	int relx = 2*n;
	int rept = 40;
	int t_step = n*6000*relx/rept;

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

}
