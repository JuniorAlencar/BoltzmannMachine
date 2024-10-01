#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <cmath>
#include <iomanip>  // Necessário para setprecision
#include <ctime>
#include <sstream>
#include <fstream>
#include <fmt/core.h>
#include "./include/nr3.h"
#include "./include/network.h"
#include "./include/forwardmethod.h"
#include "./include/InverseMethod.h"
#include "./include/LUdcmp.h"
#include "./include/write_json.h"


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
    
 }