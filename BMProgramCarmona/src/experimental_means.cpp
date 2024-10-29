#include "experimental_means.h"
using namespace std;

#include <iostream>
#include <fstream>
#include <algorithm> // provides `count`

#include "logger.h"

void read_tidy_data(const string &file_input,
                    std::vector<std::vector<int>> &M,
                    int &nrows,
                    int &ncols)
{

    auto logger = get_logger();
    ifstream data_input(file_input.c_str());
    //---------------------------------

    string first_line;
    string a;

    getline(data_input, first_line);

    // Numero de observaveis
    ncols = count(first_line.begin(), first_line.end(), ',');

    // Contar o numero de amostrar
    nrows = 0;

    while (getline(data_input, first_line))
    {
        nrows++;
    }

    nrows++;
    logger->info("read_tidy_data: read from {}", file_input); // cols
    logger->info("read_tidy_data: n_spins: {}, n_samples: {}", ncols, nrows); // cols

    // Voltar para o inicio do arquivo
    data_input.clear();
    data_input.seekg(0, ios::beg);

    // Criar matriz para receber os valores
    //  vector<vector<int>> M(rows, vector<int>(cols));
    M.resize(nrows);
    for (auto &row : M)
    {
        row.resize(ncols);
    }

    // getline(data_input, first_line);

    for (int w = 0; w < nrows; w++)
    {
        for (int p = 0; p < ncols; p++)
        {
            if (p < ncols - 1)
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

    // fechar arquivo
    data_input.close();
}

//ncols: nspins, 
//nrows: nsamples
//M: matrix with all data (nspins X nrows)
//s: first moment exp
//ss: second moment exp
//sss: third moment exp
void compute_exp_means(const vector<vector<int>> &M,
                       const int &nrows,
                       const int &ncols,
                       vector<double> &s, 
                       vector<double> &ss,
                       vector<double> &sss)
{
    
    // number combinations s_i*s_j*s_k
    int ntriplets = ncols * (ncols - 1) * (ncols - 2) / 6;
    // number combinations s_i*s_j
    int npairs = ncols * (ncols - 1) / 2;
    
    // auxiliary matrix to calculate triplet
    vector<vector<double>> M_ss(ncols, vector<double>(ncols));
    
    try{
        s.resize(ncols);
        ss.resize(npairs);
        sss.resize(ntriplets);
    }catch (const std::bad_alloc & e){
        cerr << "some problem" << endl;
        exit(1);
    }

    // First moment (magnetization)
    for (int i = 0; i < ncols; i++)
    {
        for (int row = 0; row < nrows; row++)
        {
            s[i] += M[row][i];
        }
        s[i] /= nrows;
    }


    // Second moment
    int ind = 0;
    for (int i = 0; i < ncols - 1; i++)
    {
        for (int j = i + 1; j < ncols; j++)
        {
            for (int row = 0; row < nrows; row++)
            {
                ss[ind] += M[row][i] * M[row][j];
            }
            ss[ind] /= nrows;
            
            M_ss[i][j] = ss[ind];
			M_ss[j][i] = M_ss[i][j];
            
            ++ind;
        }
    }


    // Third Moment
    ind = 0;
    for (int i = 0; i < ncols - 2; ++i)
    {
        for (int j = i + 1; j < ncols - 1; ++j)
        {
            for (int k = j + 1; k < ncols; ++k)
            {
                for(int row=0; row<nrows; ++row){
                    sss[ind] += M[row][i] * M[row][j]*M[row][k];
                }
                sss[ind] /= nrows;
                ++ind;
            }
        }
    }

  // Cij exp -> Covariance exp
  vector<double> Cij(ncols*(ncols-1)/2, 0.0);

  // Pij exp -> Pearson's correlation
  vector<double> Pij(npairs);
  
  // Calculate Pearson and Correlation exp
  ind = 0;
	
  for (int p = 0; p < ncols-1; p++)
	{
		for (int pp = p+1; pp < ncols; pp++)
		{
			Cij[ind] = ss[ind] -s[p]*s[pp];

			Pij[ind] = Cij[ind]/sqrt((1 - pow(s[p], 2))*(1 - pow(ss[pp], 2)));
			
			ind++;
		}
	}
    // Tijk = Triplet
    vector<double> Tijk(ntriplets, 0.0);

    // Triplet
    ind = 0;
    for (int i = 0; i < ncols-2; i++)
    {
        for (int j = i+1; j < ncols-1; j++)
        {
            for (int k = j+1; k < ncols; k++)
            {
                Tijk[ind] = sss[ind] - s[i]*M_ss[j][k] - s[j]*M_ss[i][k] 
                                - s[k]*M_ss[i][j] + 2*s[i]*s[j]*s[k];
                
                                
                ind++;
            }
        }
    }
}
