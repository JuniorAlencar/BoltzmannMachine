#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <cmath>
#include <iomanip>  // Necessário para setprecision
/******************************************************************************

 Codigo para tratar os dados fornecidos e transforma-los em um arquivo para
 se calcular as magnetizações e correlações.  
 
******************************************************************************/

using namespace std;

int main (int argc, char *argv[]){	
    //---------------------------------
    //Configuração para receber os da#include <iomanip>  // Necessário para std::setprecisiondos
	string file_name_input = argv[1];
	double min_erro_j	= std::stod(argv[2]);
	double min_erro_h	= std::stod(argv[3]);
	int multiply_teq 	= std::stoi(argv[4]);
	int multiply_relx 	= std::stoi(argv[5]);
	bool use_exact = (std::string(argv[6]) == "true");

    if (argc < 7) {
        std::cerr << "Uso: " << argv[0] << " <param1> <min_erro_j> <min_erro_h> <exact_solutions>" << std::endl;
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
	// Converter to string mantendo notação cientifica
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

	string file_input = "../Data/TidyData/" + file_name_input + ".dat";
    ifstream data_input (file_input.c_str());
    
	// Mag_corr
	string file_name_output;


	// nome dos arquivos para salvar----------------------
	string file_name_Cij;
	string file_name_Pij;
	string file_name_sisj;
    string file_name_Tijk;
	string file_name_sisjsk;
	string file_name_mi;
	
	if(use_exact == true){
		file_name_output = "../Data/Mag_Corr/mag_corr_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_"+ multi_teq_str + "_mrelx_" + multi_relx_str + "_exact.dat";
		file_name_Cij = "../Results/SeparateData/Cij-exp/Cij_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
		file_name_Pij = "../Results/SeparateData/Pij-exp/Pij_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
		file_name_sisj = "../Results/SeparateData/sisj-exp/sisj_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
    	file_name_Tijk = "../Results/SeparateData/Tijk-exp/Tijk_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
		file_name_sisjsk = "../Results/SeparateData/sisjsk-exp/sisjsk_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
		file_name_mi = "../Results/SeparateData/mi-exp/mi_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	}
	else{
	file_name_output = "../Data/Mag_Corr/mag_corr_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + "_metropolis.dat";
	file_name_Cij = "../Results_Metropolis/SeparateData/Cij-exp/Cij_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_Pij = "../Results_Metropolis/SeparateData/Pij-exp/Pij_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_sisj = "../Results_Metropolis/SeparateData/sisj-exp/sisj_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_Tijk = "../Results_Metropolis/SeparateData/Tijk-exp/Tijk_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_sisjsk = "../Results_Metropolis/SeparateData/sisjsk-exp/sisjsk_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_mi = "../Results_Metropolis/SeparateData/mi-exp/mi_exp_" + file_name_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	}
	
	// ---------------------
	
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

    cout << "Numero de observaveis = " << number << endl;
    cout << "Numero de amostrar = " << m << endl;

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

    //-------------------------------------------------------------------------
    //Calcular as magnetizações, segundos momentos e terceiros momentos

	int n_triplet = N*(N-1)*(N-2)/6;
	int n_duplet  = N*(N-1)/2;

    vector<double> s(N, 0.0), ss(n_duplet, 0.0), sss(n_triplet, 0.0);
	vector<vector<double>> M_ss(N, vector<double>(N));
	
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

//-----------------------------------------------------------------------------
//Salvar os dados separadamente

    //Arquivo para Cij -----------------
    

    //abrir arquivo
    ofstream file_mi (file_name_mi.c_str());

    for (int i = 0; i < N; i++)
    {
        file_mi << s[i] << endl;
    }

    //fechando arquivos
    file_mi.close();

    //-----------------------------------------------------
    //Arquivo para Cij -----------------
	//abrir arquivo
    ofstream file_Cij (file_name_Cij.c_str());

	//Arquivo para Pij------------------
	
	//Abrir arquivo
	ofstream file_Pij (file_name_Pij.c_str());

	//Arquivo para sisj
	

	//abrir arquivo
	ofstream file_sisj (file_name_sisj.c_str());

    //salvar arquivos
    for (int i = 0; i < (N*(N-1)/2); i++)
    {
        file_Cij << C[i] << endl;
		file_Pij << Pearson[i] << endl;
		file_sisj << ss[i] << endl;
    }

    //fechando arquivos
    file_Cij.close();
	file_Pij.close();
	file_sisj.close();


	//-----------------------------------------------------
    //Arquivo para Cij -----------------
    

    //abrir arquivo
    ofstream file_Tijk (file_name_Tijk.c_str());

	//Arquivo para sisjsk
	

	//abrir arquivo
	ofstream file_sisjsk (file_name_sisjsk.c_str());

    //salvar arquivos
    for (int i = 0; i < n_triplet; i++)
    {
        file_Tijk << Triplet[i] << endl;
		file_sisjsk << sss[i] << endl;
    }

    //fechando arquivos
    file_Tijk.close();
	file_sisjsk.close();


    return 0;
}