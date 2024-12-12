#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <filesystem>
#include "nr3.h"
#include "network.h"
#include "forwardmethod.h"
#include "InverseMethod.h"
#include "LUdcmp.h"
#include "create_folders.h"
#include "json_functions.h"



//#include "exp_means.h"


int main(int argc, char *argv[]){
    // set variables
 	string filename	= argv[1];
	double min_erro_j	= std::stod(argv[2]);
	double min_erro_h	= std::stod(argv[3]);
	int multiply_teq 	= std::stoi(argv[4]);
	int multiply_relx 	= std::stoi(argv[5]);
	string method = argv[6];
	
    if (argc < 7) {
        std::cerr << "Uso: " << argv[0] << " <param1> <min_erro_j> <min_erro_h> <multi_teq> <multi_relx> <exact_solutions>" << std::endl;
        return 1;
    }
    // print and check if values are available----------------------
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
    // Converter variables input to string
    ostringstream os_j;
    os_j << scientific << setprecision(2) << min_erro_j;

	ostringstream os_h;
    os_h << scientific << setprecision(2) << min_erro_h;
	
	string min_erro_j_str = os_j.str();
	string min_erro_h_str = os_h.str();
	// Create folders -------------------------------------------------------------
    c_folders cf;
    string ep;
    ep = cf.create_folders(filename, multiply_teq, multiply_relx, method, 0);
    exp_mean_calculate emc;
    
    // If file experimental means exist, open it. 
    // Else, create it and allocate in experimental means
    exp_means experimental_means;
    experimental_means = emc.exp_calculate(filename);

	// Gaussian Parameters ------------------------------------------------------
	int n;
	double mean = 1;
	double sigma = 0.25;
	double k = 10;
	//int type;
	//int H;

    // That variable are modify to save each file
    string chameleon_name_file;
    
    chameleon_name_file = "../Results/" + filename + "_corr.dat";
    
    cout << "BMfinal ..." << endl;
	
	ifstream rede (chameleon_name_file.c_str());
	
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
	
	//rede.close();
	
//-----------------------------------------------------------------------------
    //Abrir arquivo da rede
    ep = cf.create_folders(filename, multiply_teq, multiply_relx, method, 3);
    chameleon_name_file = ep + "/min_j_" + min_erro_j_str + "_min_h_" + min_erro_h_str + ".dat";
	
	
	ifstream network_in (chameleon_name_file.c_str());
	
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
	//rede para ser atualizada
	//Rede bm(n, mean, sigma, k, 0, 1);

	VecDoub bm_av_s(n, 0.0), bm_av_ss(n*(n-1)/2, 0.0);

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
	ep = cf.create_folders(filename, multiply_teq, multiply_relx, method, 2);
    chameleon_name_file = ep + "/min_j_" + min_erro_j_str + "_min_h_" + min_erro_h_str + ".dat";
	ofstream erros (chameleon_name_file.c_str());
	
	erros << "inter" << " " <<  "erroJ" << " " << "erroh" << endl; 
	
	while ((erroJ > min_erro_j || erroh > min_erro_h) && inter <= inter_max)    //(inter <= inter_max)
	{	
		srand(time(NULL)*time(NULL));

		erroJ = erroh = 0;

		eta_J = pow(inter, -0.4);
		eta_h = 2*pow(inter, -0.4);

		if (method == "exact")
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
		
		erroJ = sqrt(erroJ/bm.nbonds);
		erroh = sqrt(erroh/bm.n);

		//Salvando Erros
		erros << inter << " " << setprecision(13) << erroJ << " " << setprecision(13) << erroh << endl; 
		
		if (inter%cort == 0 || (erroJ < min_erro_j && erroh < min_erro_h))
		{
			std::cout << filename << " " << inter << " "
              << "err_J" << " " << left << setw(13) << scientific << setprecision(6) << erroJ << " "
              << "err_h" << " " << left << setw(13) << scientific << setprecision(6) << erroh << '\n';
		}			

		inter++;
	
	}
	
	//Fechar arquivos dos erros salvos 
	erros.close();
	
//-----------------------------------------------------------------------------
//Arquivo para salvar a rede obtida	

    //Abrir arquivo da rede
    ep = cf.create_folders(filename, multiply_teq, multiply_relx, method, 3);
    chameleon_name_file = ep + "/min_j_" + min_erro_j_str + "_min_h_" + min_erro_h_str + ".dat";

	ofstream network (chameleon_name_file.c_str());

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
	
    chameleon_name_file = "../Results/" + filename + "_corr.dat";
    ofstream mag_corr (chameleon_name_file.c_str());

	VecDoub bm_C(n*(n-1)/2, 0.0);

	//correlação pearson ()
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
    
	//-----------------------------------------------------
	//Calcular os tripletos

	int n_triplet = n*(n-1)*(n-2)/6;
	vector<double> bm_av_sss(n_triplet), ising_Triplet(n_triplet);
	vector<vector<double>> ising_M_av_ss(n, vector<double>(n));

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
				ind++;
			}
		}
	}
	return 0;
}

