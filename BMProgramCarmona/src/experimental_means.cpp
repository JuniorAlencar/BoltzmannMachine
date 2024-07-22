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

void compute_exp_means(const vector<vector<int>> &M,
                       const int &nrows,
                       const int &ncols,
                       vector<double> &s, vector<double> &ss,
                       vector<double> &sss)
{

    int ntriplets = ncols * (ncols - 1) * (ncols - 2) / 6;
    int npairs = ncols * (ncols - 1) / 2;
	
    try{
        s.resize(ncols);
        ss.resize(npairs);
        sss.resize(ntriplets);
    }catch (const std::bad_alloc & e){
        cerr << "some problem" << endl;
        exit(1);
    }

    // Primeiro momento (magnetização)
    for (int i = 0; i < ncols; i++)
    {
        for (int row = 0; row < nrows; row++)
        {
            s[i] += M[row][i];
        }
        s[i] /= nrows;
    }



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
            ++ind;
        }
    }



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


}
