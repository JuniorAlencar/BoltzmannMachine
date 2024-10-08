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
#include "nr3.h"
#include "network.h"
#include "LUdcmp.h"
#include "forwardmethod.h"
#include "InverseMethod.h"
#include "json_functions.hpp"
#include "exp_means.hpp"


int main(int argc, char *argv[]){
    // set variables
 	string filename	= argv[1];
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
    
    exp_mean_calculate exp_funcs;
    // If file experimental means exist, open it. 
    // Else, create it and allocate in experimental means
    exp_means experimental_means;
    experimental_means = exp_funcs.exp_calculate(filename);
    
    return 0;
}

