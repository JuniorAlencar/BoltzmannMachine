#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <fstream>
#include <fmt/core.h>
#include "./include/nr3.h"
#include "./include/network.h"
#include "./include/LUdcmp.h"
#include "./include/forwardmethod.h"
#include "./include/InverseMethod.h"
#include "./include/write_json.h"
#include "./include/read_json.h"
#include <filesystem>

using namespace std;
using json = nlohmann::json;


 int main(int argc, char *argv[]){
    // set variables
 	string text_name	= argv[1];
	double min_erro_j	= std::stod(argv[2]);
	double min_erro_h	= std::stod(argv[3]);
	int multiply_teq 	= std::stoi(argv[4]);
	int multiply_relx 	= std::stoi(argv[5]);
	bool use_exact = (std::string(argv[6]) == "true");
	
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
    
    // Count the number of spins in sample ---------------------------------
    string file_input = "../Data/TidyData/" + text_name + ".dat";
    ifstream data_input (file_input.c_str());
    string first_line;
    string a;

    getline(data_input, first_line);
    int N = count(first_line.begin(), first_line.end(), ',');
    
    // Number of spins
    int nspins = N + 1; 
    
    // Converter variables to string ------------------------------------------	
	ostringstream os_teq;
    // teq = nspins * multiply_teq
	os_teq << multiply_teq*nspins;

	ostringstream os_relx;
    // relx = nspins * multiply_relx
	os_relx << multiply_relx*nspins;
	
	string teq_str = os_teq.str();
	string relx_str = os_relx.str();
	
	ostringstream os_j;
    // Converver jmin to string keeping scientific notation
    os_j << scientific << setprecision(2) << min_erro_j;

	ostringstream os_h;
    // Converver hmin to string keeping scientific notation
    os_h << scientific << setprecision(2) << min_erro_h;

	string min_erro_j_str = os_j.str();
	string min_erro_h_str = os_h.str();
	
    // count the number of samples
    int m = 0;
	
	int number = N + 1;
    
	while (getline(data_input, first_line))
    {
        m++;
    }

	m++;

    cout << "Number of spins = " << number << endl;
    cout << "Number of samples = " << m << endl;

    // return to the beginning of the file
    data_input.clear();
	data_input.seekg (0, ios::beg);

    
	// Create the matrix to receive the values
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
            // Update matrix
			M[w][p] = stoi(a);				
        }

	}

    // Close file
    data_input.close();
    
    //===>>> CALCULATE ALL PROPERTIES EXPERIMENTAL <<<===

	
	// Check if file exp exist
	string file_means = use_exact ? "../Results/" + text_name + "_exp.json" : "../Results_Metropolis/" + text_name + "_exp.json";
	
	// load struct with means exp values
	exp_means my_means;
	// json struct with means to append with ising means
	json json_means_exp;
	
	// If exist, just open it
	if (filesystem::exists(file_means)) {
        cout << "file exist, opening..." << endl;
		// 
		pair<exp_means, json> result = load_json_exp(file_means);
		
		// Acessando a struct exp_means
    	my_means = result.first;
    
    	// Acessando o objeto JSON
    	json_means_exp = result.second;
    } 
	
	// Else, calculate and create it
	else {
        cout << "file means don't exist, calculate..." << endl;
		// First, second and third moment
		vector<double> s(N, 0.0), ss(n_duplet, 0.0), sss(n_triplet, 0.0);  
		// Auxiliary matrix to calculate the triplet and third moment
		vector<vector<double>> M_ss(N, vector<double>(N));
		
		// First moment (magnetization)
		for (int p = 0; p < N; p++)
		{
			for(int w = 0; w < m; w++)
			{
				s[p] += M[w][p];
			}
			
			s[p] /= m;
		}

		// Second moment
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
		// Third moment
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

		// Triplet
		vector<double> Triplet(n_triplet, 0.0);
		
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

		// Covariance experimental (Cij)
		vector<double> C(N*(N-1)/2, 0.0);
		
		// Correlation experimental (Pij)
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
		
		// Update struct exp_means with values
		my_means.av_s = s;
		my_means.av_ss = ss;
		my_means.av_sss = sss;
		my_means.Cij_exp = C;
		my_means.Pij_exp = Pearson;
		my_means.Tijk_exp = Triplet;
		
		// Json saved and created
		json_means_exp = create_json_exp(my_means, text_name);
		
		cout << "json " << endl;
    }
	
    
    
    
	// //===>>> CALCULATE ALL PROPERTIES ISING <<<===
	// cout << "start ising means..." << endl;

	// // Gaussian Parameters
	// int n;
	// double mean = 1;
	// double sigma = 0.25;
	// double k = 10;
	// int type;
	// int H;
	
	// string file_rede_input = use_exact ? "exact" : "metropolis";

	// ifstream rede (file_rede_input.c_str());
	
	// rede >> n;
	
	// VecDoub av_s(n, 0.0), av_ss(n*(n-1)/2, 0.0), C(n*(n-1)/2, 0.0);
	
	// //Rede r(tamanho, media, desvio, k, type = 0(random) 1(tree), H = 0(no field) 1(ConstField) -1(ConstField) 2(RandField))
	// Rede r(n, mean, sigma, k, 0, 1);	
	
	// for (int i = 0; i < r.nbonds; i++)
	// {
	// 	rede >> av_ss[i] >> C[i];
		
	// 	if (i < r.n)
	// 		rede >> av_s[i];	
		 
	// }


 }
