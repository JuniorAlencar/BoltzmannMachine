#ifndef EXP_MEANS_H
#define EXP_MEANS_H

#include <iostream>
#include <ctime>
#include <math.h>
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <fstream>
#include <boost/filesystem.hpp>

//#include "./create_folders.h"


using namespace std;
// nspins: Number of spins
// t_eq: Multiply_teq*nspins
// relx: Multiply_relx*nspins
// bm_av_s: First moment experimental
// bm_av_ss: Second moment experimental
// bm_av_sss: Third moment experimental
// C_exp: Covariance experimental
// Pij: Correlation experimental
// Tijk: Triplet experimental
struct exp_means{
    int nspins;
    int t_eq;
    int relx;
    vector<double> bm_av_s;           
    vector<double> bm_av_ss;          
    vector<double> bm_av_sss;      
    vector<double> C_exp;      
    vector<double> Pij_exp;    
    vector<double> Tijk_exp;
};


exp_means experimental_means(const string &filename, 
                             const int &multiply_teq, 
                             const int &multiply_relx,
                             const double &min_err_h,
                             const double &min_err_j,
                             const bool &use_exact
                        ){

    struct exp_means MV;

    string file_input = "../Data/TidyData/" + filename + ".dat";
    ifstream data_input (file_input.c_str());


    string first_line;
    string a;

    getline(data_input, first_line);

    //Numero de observaveis
    int N = count(first_line.begin(), first_line.end(), ',');

    //Contar o numero de amostrar
    int m = 0;

    int number = N + 1;

    while (getline(data_input, first_line))
    {
        m++;
    }

    m++;

    //cout << "Numero de observaveis = " << number << endl;
    //cout << "Numero de amostrar = " << m << endl;

    //Voltar para o inicio do arquivo
    data_input.clear();
    data_input.seekg (0, ios::beg);


    //Criar matriz para receber os valores
    vector<vector<int>> M(m, vector<int>(N));


    //getline(data_input, first_line);

    for (int w = 0; w < m; w++)
    {
        for(int p = 0; p < N; p++)
        {
            if (p < N - 1)
            {
                getline(data_input, a, ',');
            }
            else
            {
                getline(data_input, a);
            }           

            M[w][p] = stoi(a);				
        }

    }

    //fechar arquivo
    data_input.close();

    int real_teq = multiply_teq * number;
    int real_relx = multiply_relx * number;
    int n_triplet = N*(N-1)*(N-2)/6;
    int n_duplet  = N*(N-1)/2;

    vector<double> s(N, 0.0), ss(n_duplet, 0.0), sss(n_triplet, 0.0);
    vector<vector<double>> M_ss(N, vector<double>(N));

    //magnetizações, segundos momentos e terceiros momentos

    //Primeito momento (magnetização)
    for (int p = 0; p < N; p++)
    {
        for(int w = 0; w < m; w++)
        {
            s[p] += M[w][p];
        }
        
        s[p] /= m;
    }

    //Segundo momento
    int ind = 0;
    for (int p = 0; p < N-1; p++)
    {
        for(int pp = p+1; pp < N; pp++)
        {		
            for (int w = 0; w < m; w++)
            {
                ss[ind] += M[w][p]*M[w][pp];
            }
        
            ss[ind] /= m;

            M_ss[p][pp] = ss[ind];
            M_ss[pp][p] = M_ss[p][pp];
        
            ind++;
        }
    }

    //Covariancia
    vector<double> C(N*(N-1)/2, 0.0);

    // Corelação de Pearson
    vector<double> Pearson(n_duplet);

    ind = 0;
    for (int p = 0; p < N-1; p++)
    {
        for (int pp = p+1; pp < N; pp++)
        {
            C[ind] = ss[ind] - s[p]*s[pp];

            Pearson[ind] = C[ind]/sqrt((1 - pow(s[p], 2))*(1 - pow(s[pp], 2)));
            
            ind++;
        }
    }

    //Terceiro momento
    ind = 0;
    for (int i = 0; i < N-2; i++)
    {
        for (int j = i+1; j < N-1; j++)
        {
            for (int k = j+1; k < N; k++)
            {
                for (int s = 0; s < m; s++)
                {
                    sss[ind] += M[s][i]*M[s][j]*M[s][k];
                }
                
                sss[ind] /= m;
                ind++;
            }
        }
    }

    vector<double> Triplet(n_triplet, 0.0);

    //Terceiro momento centrado (tripleto)
    ind = 0;
    for (int i = 0; i < N-2; i++)
    {
        for (int j = i+1; j < N-1; j++)
        {
            for (int k = j+1; k < N; k++)
            {
                Triplet[ind] = sss[ind] - s[i]*M_ss[j][k] - s[j]*M_ss[i][k] 
                                - s[k]*M_ss[i][j] + 2*s[i]*s[j]*s[k];
                
                                
                ind++;
            }
        }
    }
    
    MV.nspins = number;    // number of spins
    MV.t_eq = real_teq;    // t_eq value
    MV.relx = real_relx;   // relx value
    MV.bm_av_s = s;        // First moment exp
    MV.bm_av_ss = ss;      // Second moment exp
    MV.bm_av_sss = sss;    // Third moment exp
    MV.C_exp = C;          // Covariance exp
    MV.Pij_exp = Pearson;  // Correlation exp
    MV.Tijk_exp = Triplet; // Triplet exp

    string file_name_output;	
    // converter t_eq and relx to string=================
    ostringstream os_teq;
	os_teq << real_teq;

	ostringstream os_relx;
	os_relx << real_relx;

    string teq_str = os_teq.str();
	string relx_str = os_relx.str();

	// Converter min_err_j and min_err_h to string=============
	ostringstream os_j;
    os_j << std::scientific << setprecision(2) << min_err_j;

	ostringstream os_h;
    os_h << std::scientific << setprecision(2) << min_err_h;

	string min_erro_j_str = os_j.str();
	string min_erro_h_str = os_h.str();
    
    if(use_exact==true)
		file_name_output = "../Data/Mag_Corr/" + filename + "/exact/teq_" + teq_str + 
		"/relx_" + relx_str + "/exp_j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + ".dat";
	
	else
		file_name_output = "../Data/Mag_Corr/" + filename + "/metropolis/teq_" + teq_str + 
		"/relx_" + relx_str + "/exp_j_min_" + min_erro_j_str + "_h_min_" + min_erro_h_str + ".dat";

    //Cria arquivo para salvar correlação e magnetização-----------------------
	ofstream mag_corr (file_name_output.c_str());
	
	mag_corr << N << endl;
	
	for (int i = 0; i < N*(N-1)/2; i++)
	{
		mag_corr << " " << ss[i] << " " << C[i];
		
		if (i < N)
			mag_corr << " " << s[i];
			
		mag_corr << endl;
	}

    mag_corr.close();
    
    return MV;
}


#endif /* _EXP_MEANS_H_ */