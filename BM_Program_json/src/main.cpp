#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <fstream>
#include "nr3.h"
#include "network.h"
#include "LUdcmp.h"
#include "forwardmethod.h"
#include "InverseMethod.h"
#include "exp_means.h"


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
    c_folders cf;
    int tp = 0;
    string ep;
    ep = cf.create_folders(filename, multiply_teq, multiply_relx, method, tp);

    exp_mean_calculate emc;
    // If file experimental means exist, open it. 
    // Else, create it and allocate in experimental means
    exp_means experimental_means;
    experimental_means = emc.exp_calculate(filename);
    

	// Gaussian Parameters
	int n;
	double mean = 1;
	double sigma = 0.25;
	double k = 10;
	int type;
	int H;

    cout << "BMfinal ..." << endl;
	
    string file_rede_input = "../Results/" + filename + "_mag_corr.dat";
	
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
		 
	};
    
    return 0;
}

