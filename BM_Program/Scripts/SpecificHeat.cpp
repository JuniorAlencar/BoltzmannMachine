#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include "./include/nr3.h"
#include <cmath>
#include "./include/network.h"
#include "./include/forwardmethod.h"

using namespace std;

int main(int argc, char *argv[]){


	string text_input = argv[1];
	double min_erro_j	= std::stod(argv[2]);
	double min_erro_h	= std::stod(argv[3]);
	bool use_exact = (std::string(argv[4]) == "true");

    if (argc < 5) {
        std::cerr << "Uso: " << argv[0] << " <param1> <min_erro_j> <min_erro_h> <exact_solutions>" << std::endl;
        return 1;
    }

    try {
        double min_erro_j = std::stod(argv[2]);
        double min_erro_h = std::stod(argv[3]);

        std::cout << "min_erro_j: " << min_erro_j << std::endl;
        std::cout << "min_erro_h: " << min_erro_h << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Argumento inválido: " << e.what() << std::endl;
        return 1;
    } catch (const std::out_of_range& e) {
        std::cerr << "Valor fora do intervalo: " << e.what() << std::endl;
        return 1;
    }


	int n;

	int cont;
	int cort = 100;
	// Convertendo os erros para string
	std::ostringstream os_j;
    os_j << std::scientific << std::setprecision(2) << min_erro_j;

	std::ostringstream os_h;
    os_h << std::scientific << std::setprecision(2) << min_erro_h;

	std::string min_erro_j_str = os_j.str();
	std::string min_erro_h_str = os_h.str();


//-----------------------------------------------------------------------------
//Abrir arquivo da rede

	string file_network_name;
	string file_cap_linear_name;
	string file_energia_name;
	string file_mag_vs_T_name;
	string file_susc_vs_T_name;
// nome dos arquivos-------------------------------------------------
	
	if(use_exact == true){
		// Network file to exact solutions
		file_network_name = "../Results/Network/network_" + text_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + ".dat";
		file_cap_linear_name = "../Results/SpecificHeat/linear_specific_heat_" + text_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + ".dat";
		file_energia_name = "../Results/Energy/energy_" + text_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + ".dat";
		file_mag_vs_T_name = "../Results/Magnetization_vs_T/mag_T_" + text_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + ".dat";
		file_susc_vs_T_name = "../Results/Magnetization_vs_T/susc_T_" + text_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + ".dat";
	}
	
	else{
	// Network file to metropolis algorithm
	file_network_name = "../Results_Metropolis/Network/network_" + text_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + ".dat";
	file_cap_linear_name = "../Results_Metropolis/SpecificHeat/linear_specific_heat_" + text_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + ".dat";
	file_energia_name = "../Results_Metropolis/Energy/energia_" + text_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + ".dat";
	file_mag_vs_T_name = "../Results_Metropolis/Magnetization_vs_T/mag_T_" + text_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + ".dat";
	file_susc_vs_T_name = "../Results_Metropolis/Magnetization_vs_T/susc_T_" + text_input + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + ".dat";
	}
	

	ifstream network (file_network_name.c_str());
	int number = n+1;
	
	network >> number;
	
	Rede r(n, 0, 0, 0, 0, 0);
	
	for(int i = 0; i < r.nbonds; i++)
	{
		network >> r.J[i];
				
		if (i < n)
			network >> r.h[i];
	
	}
	
	network.close();

//Escala linear

	//variaveis para MC
	int t_eq = n*150;
	int relx = 2*n;
	int rept = 40;
	int t_step = n*6000*relx/rept;

	//bool use_exact = false;
	
	ofstream cap_li (file_cap_linear_name.c_str());


	ofstream ene (file_energia_name.c_str());
	
	double E, E2, sp;

	//Arquivo para salvar magnetização versus temperatura
    
    
    ofstream magT (file_mag_vs_T_name.c_str());

    //Arquivo para salvar Susceptibility versus temperatura
    
    
    ofstream susc (file_susc_vs_T_name.c_str());

    double M, M2;
    double S;

	cout << "Calor Especifico dados " << text_input << endl;
	
	for(double T = 3; T >= 0.01; T -= 0.01)
	{
		E = E2 = 0;
		M = M2 =0;

		metropolis_cap_mag (r, t_eq, t_step, relx, rept, T, E, E2, M, M2);
	
		
		//metropolis_cap (r, t_eq, t_step, relx, rept, T, E, E2);
		
		sp = pow(1.0/T, 2)*(E2 - pow(E, 2))/r.n;
		
		cap_li << T << " " << sp << endl;
		ene << T << " " << E << endl;

		S = n*(M2 - M*M)/T;

		magT << T << " " << abs(M)/n << endl;
        susc << T << " " << S << endl;

		//cout << text_input << " " << left << setw(10) << T << sp << endl;
	}
	
	cap_li.close();
	ene.close();
	magT.close();
	susc.close();
	
//-----------------------------------------------------------------------------
/*
//-----------------------------------------------------------------------------
//Escala logaritmica

	string file_cap_log_name = "./results/log_specific_heat/log_specific_heat_" + text_input;
	
	ofstream cap_log (file_cap_log_name.c_str());
	
	//double E, E2, sp;
	double T;
	
	for(double x = -1; x <= 2; x += 0.01)
	{
	
		T = pow(10, x);		
			
		E = E2 = 0;
	
		exact_solution_cap (r, E, E2, 1.0/T);

		//metropolis_cap (r, t_eq, t_step, relx, rept, T, E, E2);
		
		sp = pow(1.0/T, 2)*(E2 - pow(E, 2))/r.n;
		
		cap_log << left << setw(6) << T << " " << left << setw(10) << sp << endl;
		
		//cout << left << setw(10) << T << sp << endl;
	}
	
	cap_log.close();

*/

	
//-----------------------------------------------------------------------------








	return 0;
}
