#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <fmt/core.h>
#include <chrono>
#include <omp.h>
#include "./include/nr3.h"
#include "./include/network.h"
#include "./include/forwardmethod.h"
#include "./include/InverseMethod.h"
#include "./include/LUdcmp.h"


using namespace std::chrono;
using namespace std;

int main(int argc, char *argv[]){
	// Gaussian Parameters
	int n;
	double mean = 0.0;
	double sigma = 0.25;
	double k = 10;
	int type = 0;
	int H = 2;
	
	string text_name	= argv[1];
	double min_erro_j	= stod(argv[2]);
	double min_erro_h	= stod(argv[3]);
	int multiply_teq 	= stoi(argv[4]);
	int multiply_relx 	= stoi(argv[5]);
	int seed 			= stoi(argv[6]);
	string method = argv[7];
	
    if (argc < 8) {
        cerr << "Uso: " << argv[0] << " <filename> <min_erro_j> <min_erro_h> <multi_teq> <multi_relx> <seed> <method>" << endl;
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
	
	// Converter para string
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
	
	// Filenames to results
	string file_name_erros;
	string file_mag_corr_output;
	string file_name_Jij;
	string file_name_Cij;
	string file_name_Pij;
	string file_name_sisj;
	string file_name_Tijk;
	string file_name_sisjsk;
	string file_name_CorrJij;
	string file_name_PJij;
	string file_name_hi;
	string file_name_mi;

	string file_rede_input = "../Data/Mag_Corr/mag_corr_exp_" + text_name + ".dat";
	string file_network_name = "../Results/" + method +  "/Network/network_" + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	string file_rede_output = "../Results/" + method +  "/Network/network_" + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	
	// files name
	file_name_erros = "../Results/" + method +  "/Erro/erro_" + text_name  + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + "_seed_" + to_string(seed) + ".dat";
	file_mag_corr_output = "../Results/" + method +  "/Mag_Corr_ising/mag_corr_ising_" + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_Jij = "../Results/" + method +  "/SeparateData/Jij/Jij_" +  text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_Cij = "../Results/" + method +  "/SeparateData/Cij-ising/Cij_ising_" + text_name  + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_Pij = "../Results/" + method +  "/SeparateData/Pij-ising/Pij_ising_"  + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_sisj = "../Results/" + method +  "/SeparateData/sisj-ising/sisj_ising_" + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_Tijk = "../Results/" + method +  "/SeparateData/Tijk-ising/Tijk_ising_" + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_sisjsk = "../Results/" + method +  "/SeparateData/sisjsk-ising/sisjsk_ising_" + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_CorrJij = "../Results/" + method +  "/CorrJij/CorrJij_" + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_PJij = "../Results/" + method +  "/PJij/PJij_" + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_hi = "../Results/" + method +  "/SeparateData/hi/hi_" + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	file_name_mi = "../Results/" + method +  "/SeparateData/mi-ising/mi_ising_" + text_name + "_err_j_" + min_erro_j_str + "_err_h_" + min_erro_h_str + "_mteq_" + multi_teq_str + "_mrelx_" + multi_relx_str + ".dat";
	string file_complex_time = "../Results/" + method + "/Complex_time/" + "time_method_" + to_string(seed);
	//-----------------------------------------------------------------------------
	//Ler o arquivo com as correlações e magnetizações de um certo arquivo

	//create_folders();
	
	cout << "BMfinal ..." << endl;
	
	ifstream rede (file_rede_input.c_str());
	
	rede >> n;
	// Generator to MC
	std::mt19937 gen(seed); 
	VecDoub av_s(n, 0.0), av_ss(n*(n-1)/2, 0.0), C(n*(n-1)/2, 0.0);
	
	//Rede r(tamanho, media, desvio, k, type = 0(random) 1(tree), H = 0(no field) 1(ConstField) -1(ConstField) 2(RandField))
	Rede r(n, mean, sigma, k, 0, H, gen);	
	
	for (int i = 0; i < r.nbonds; i++)
	{
		rede >> av_ss[i] >> C[i];
		
		if (i < r.n)
			rede >> av_s[i];	
		 
	}
	
	
	
	//rede.close();
	
//-----------------------------------------------------------------------------
//Abrir arquivo da rede

	ifstream network_in (file_network_name.c_str());
	
	network_in >> r.n;
	
	//rede para ser atualizada
	Rede bm(n, 0, 0, 0, 0, 0, gen);
	
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

	VecDoub_IO bm_av_s(n, 0.0), bm_av_ss(n*(n-1)/2, 0.0);

	//variaveis para MC
	int t_eq = n*multiply_teq; // 150
	int relx = n*multiply_relx; // 2
	int rept = 60;
	//int t_step = n*6000*relx/rept;
	int t_step = n*6000*relx/rept;
	
	double erroJ = 1, erroh = 1;
	double dJ, dh;
	int cort = 1000;
	int inter = 1;//= inter_ini;
	int inter_max = 300000;
	
	double eta_J = 0.05;//atof(argv[2]);
	double eta_h = 0.03;

	//Arquivo para salvar os erros ao longo do tempo


	ofstream erros (file_name_erros.c_str());

	erros << "inter" << " " <<  "erroJ" << " " << "erroh" << endl; 
	
	// const double J_MAX = 3.0;   // Valor limite para J
	// const double H_MAX = 3.0;   // Valor limite para h
	
	std::ofstream times(file_complex_time);
	if(method == "parallel_tempering")
		times << "time(s) " << "swap_accept_rate" << endl;
	else
		times << "time(s)" << endl;
	// Running with inter_max interations
	// while (inter <= inter_max)
	// {	

	// 	erroJ = erroh = 0;

	// 	eta_J = pow(inter, -0.4);
	// 	eta_h = 2*pow(inter, -0.4);

	// 	if(method == "metropolis"){
	// 		auto start = high_resolution_clock::now();
	// 		metropolis_bm(bm, bm_av_s, bm_av_ss, t_eq, t_step, relx, rept, 1, gen);
	// 		auto end = high_resolution_clock::now();
	// 		std::chrono::duration<double> duration_sec = end - start;
	// 		times << duration_sec.count() << std::endl;
	// 	}
		
	// 	else if(method == "exact" && n < 25){
	// 		auto start = high_resolution_clock::now();
	// 		exact_solution_bm (bm, bm_av_s, bm_av_ss, 1);
	// 		auto end = high_resolution_clock::now();
	// 		std::chrono::duration<double> duration_sec = end - start;
	// 		times << duration_sec.count() << std::endl;
	// 	}
		
	// // else if (method == "parallel_tempering") {
	// // 	// eta_J = pow(inter, -0.5);
	// // 	// eta_h = 2*pow(inter, -0.5);
	// // 	// --- Cooling + Reheating cycles ---
	// // 	int cycle_length = 100;  // Número de iterações por ciclo
	// // 	double cooling_factor = pow(inter, -0.5);  // Decaimento normal

	// // 	// Reheating multiplier: aumenta temporariamente o passo em ciclos ímpares
	// // 	double reheating_multiplier = 1.0;
	// // 	if ( (inter / cycle_length) % 2 == 1) {
	// // 		reheating_multiplier = 2.0;  // Pode usar 1.5 se quiser mais suave
	// // 	}

	// // 	eta_J = cooling_factor * reheating_multiplier;
	// // 	eta_h = 2 * cooling_factor * reheating_multiplier;

	
	// // 	double T_min, T_max;
	// // 	std::vector<double> temperatures, energy_per_replica;
	// // 	double swap_acceptance_ratio = 0.0;
	// // 	int n_replicas = 12;
	
	// // 	// Se for a primeira interação: temperatura padrão
	// // 	if (inter == 1) {
	// // 		T_min = 0.2;
	// // 		T_max = T_min * 20.0;  // ou outro multiplicador, se quiser
	// // 		temperatures.resize(n_replicas);
	// // 		for (int i = 0; i < n_replicas; ++i)
	// // 			temperatures[i] = T_min * pow(T_max / T_min, (double)i / (n_replicas - 1));
	// // 	}
	// // 	else {
	// // 		// Após a primeira interação: adaptar com energy_per_replica
	// // 		double mean_E = 0.0, var_E = 0.0;
	// // 		for (double E : energy_per_replica) mean_E += E;
	// // 		mean_E /= n_replicas;
	// // 		for (double E : energy_per_replica) var_E += (E - mean_E) * (E - mean_E);
	// // 		var_E /= n_replicas;
	// // 		double sigma_E = sqrt(var_E);
	
	// // 		if (sigma_E > 1e-8)
	// // 			T_min = 1.0 / (2.0 * sigma_E);
	// // 		else
	// // 			T_min = 2.0;  // fallback para caso sigma_E muito pequeno
	
	// // 		T_max = T_min * 10.0;  // ou 7.0, 20.0, conforme preferir
	
	// // 		temperatures.resize(n_replicas);
	// // 		for (int i = 0; i < n_replicas; ++i)
	// // 			temperatures[i] = T_min * pow(T_max / T_min, (double)i / (n_replicas - 1));
	// // 	}
	
	// // 	auto start = high_resolution_clock::now();
	
	// // 	parallel_tempering_multi(
	// // 		bm, n_replicas, T_min, T_max,
	// // 		bm_av_s, bm_av_ss,
	// // 		t_eq, t_step, relx, rept,
	// // 		n, mean, sigma,
	// // 		type, H, energy_per_replica, temperatures,
	// // 		swap_acceptance_ratio, gen
	// // 	);
	
	// // 	auto end = high_resolution_clock::now();
	// // 	std::chrono::duration<double> duration_sec = end - start;
	// // 	times << duration_sec.count() << " " << swap_acceptance_ratio << std::endl;
	
	// // 	// Diagnóstico no terminal
	// // 	std::cout << "Swap acceptance rate: " << swap_acceptance_ratio << std::endl;
	// // 	for (int i = 0; i < n_replicas; ++i)
	// // 		std::cout << "T = " << temperatures[i] << ", <E> = " << energy_per_replica[i] << std::endl;
	// // }
	// else if (method == "parallel_tempering") {
	// 	eta_J = pow(inter, -0.5);
	// 	eta_h = 2*pow(inter, -0.5);
	// 	// --- Cooling + Reheating cycles ---
	// 	// int cycle_length = 100;  // Número de iterações por ciclo
	// 	// double cooling_factor = pow(inter, -0.5);  // Decaimento normal
	
	// 	// // Reheating multiplier: aumenta temporariamente o passo em ciclos ímpares
	// 	// double reheating_multiplier = 1.0;
	// 	// if ( (inter / cycle_length) % 2 == 1) {
	// 	// 	reheating_multiplier = 2.0;  // Pode usar 1.5 se quiser mais suave
	// 	// }
	
	// 	// eta_J = cooling_factor * reheating_multiplier;
	// 	// eta_h = 2 * cooling_factor * reheating_multiplier;
	
	// 	// --- Temperaturas e réplicas ---
	// 	double T_min, T_max;
	// 	int n_replicas = 12;
	
	// 	// As temperaturas e as energias *devem ser mantidas entre as iterações*!
	// 	static std::vector<double> temperatures;       // Static: persiste entre as interações
	// 	static std::vector<double> energy_per_replica; // Static: idem
	
	// 	double swap_acceptance_ratio = 0.0;
	
	// 	if (inter == 1) {
	// 		T_min = 0.2;
	// 		T_max = T_min * 20.0;
	
	// 		temperatures.resize(n_replicas);
	// 		for (int i = 0; i < n_replicas; ++i)
	// 			temperatures[i] = T_min * pow(T_max / T_min, (double)i / (n_replicas - 1));
	
	// 		// Zerar energy_per_replica só no começo
	// 		energy_per_replica.assign(n_replicas, 0.0);
	// 	} 
	// 	else {
	// 		// --- Atualizar temperaturas com base nas energias ---
	// 		double mean_E = 0.0, var_E = 0.0;
	// 		for (double E : energy_per_replica) mean_E += E;
	// 		mean_E /= n_replicas;
	
	// 		for (double E : energy_per_replica)
	// 			var_E += (E - mean_E) * (E - mean_E);
	// 		var_E /= n_replicas;
	
	// 		double sigma_E = sqrt(var_E);
	
	// 		if (sigma_E > 1e-8)
	// 			T_min = 1.0 / (2.0 * sigma_E);
	// 		else
	// 			T_min = 2.0;
	
	// 		T_max = T_min * 10.0;
	
	// 		temperatures.resize(n_replicas);
	// 		for (int i = 0; i < n_replicas; ++i)
	// 			temperatures[i] = T_min * pow(T_max / T_min, (double)i / (n_replicas - 1));
	// 	}
	
	// 	// --- Executar Parallel Tempering ---
	// 	auto start = high_resolution_clock::now();
	
	// 	parallel_tempering_multi(
	// 		bm, n_replicas, T_min, T_max,
	// 		bm_av_s, bm_av_ss,
	// 		t_eq, t_step, relx, rept,
	// 		n, mean, sigma,
	// 		type, H, energy_per_replica, temperatures,
	// 		swap_acceptance_ratio, gen
	// 	);
	
	// 	auto end = high_resolution_clock::now();
	// 	std::chrono::duration<double> duration_sec = end - start;
	// 	times << duration_sec.count() << " " << swap_acceptance_ratio << std::endl;
	
	// 	// --- Diagnóstico ---
	// 	std::cout << "Swap acceptance rate: " << swap_acceptance_ratio << std::endl;
	// 	for (int i = 0; i < n_replicas; ++i)
	// 		std::cout << "T = " << temperatures[i] << ", <E> = " << energy_per_replica[i] << std::endl;
	
	// }
	
	
		
	// 	for (int i = 0; i < bm.nbonds; i++)
	// 	{
	// 		if (i < bm.n)
	// 		{
	// 			dh = eta_h*(bm_av_s[i] - av_s[i]);
	// 			erroh += pow(bm_av_s[i] - av_s[i], 2);
	// 			bm.h[i] -= dh;
	// 		}
			
	// 		dJ = eta_J*(bm_av_ss[i] - av_ss[i]);
	// 		erroJ += pow(bm_av_ss[i] - av_ss[i], 2);
	// 		bm.J[i] -= dJ;

	// 	}

		
	// 	erroJ = sqrt(erroJ/bm.nbonds);
	// 	erroh = sqrt(erroh/bm.n);

		
	// 	//Salvando Erros
	// 	erros << inter << " " << setprecision(13) << erroJ << " " << setprecision(13) << erroh << endl; 
		
	// 	if (inter%cort == 0)
	// 	{
	// 		std::cout << text_name << " " << inter << " "
    //           << "err_J" << " " << left << setw(13) << scientific << setprecision(6) << erroJ << " "
    //           << "err_h" << " " << left << setw(13) << scientific << setprecision(6) << erroh << '\n';
	// 	}			
		
	// 	// If J and h converge, break loop
	// 	if(erroJ <= min_erro_j && erroh <= min_erro_h){
	// 		cout << "Converged with: " << inter << " interations\n";
	// 		break;
	// 	}
		
	// 	inter++;
	
	// 	}
	while (inter <= inter_max) {
		erroJ = erroh = 0.0;
	
		// --- Cooling + Reheating (opcional) ---
		// int cycle_length = 100;
		// double cooling_factor = pow(inter, -0.5);
		// double reheating_multiplier = ((inter / cycle_length) % 2 == 1) ? 2.0 : 1.0;
	
		// eta_J = cooling_factor * reheating_multiplier;
		// eta_h = 2.0 * eta_J;
		eta_J = pow(inter, -0.4);
		eta_h = 2.0 * eta_J;
		if (method == "metropolis") {
			auto start = high_resolution_clock::now();
			metropolis_bm(bm, bm_av_s, bm_av_ss, t_eq, t_step, relx, rept, 1, gen);
			auto end = high_resolution_clock::now();
			std::chrono::duration<double> duration_sec = end - start;
			times << duration_sec.count() << std::endl;
		}
		else if (method == "exact" && n < 25) {
			auto start = high_resolution_clock::now();
			exact_solution_bm(bm, bm_av_s, bm_av_ss, 1);
			auto end = high_resolution_clock::now();
			std::chrono::duration<double> duration_sec = end - start;
			times << duration_sec.count() << std::endl;
		}
		else if (method == "parallel_tempering") {
			eta_J = pow(inter, -0.5);
			eta_h = 2.0 * eta_J;
			static std::vector<double> temperatures;
			static std::vector<double> energy_per_replica;
	
			double T_min, T_max;
			double swap_acceptance_ratio = 0.0;
			int n_replicas = 12;
	
			if (inter == 1) {
				T_min = 0.2;
				T_max = T_min * 20.0;
				temperatures.resize(n_replicas);
				for (int i = 0; i < n_replicas; ++i)
					temperatures[i] = T_min * pow(T_max / T_min, (double)i / (n_replicas - 1));
				energy_per_replica.assign(n_replicas, 0.0);
			} 
			else {
				double mean_E = 0.0, var_E = 0.0;
				for (double E : energy_per_replica) mean_E += E;
				mean_E /= n_replicas;
				for (double E : energy_per_replica)
					var_E += (E - mean_E) * (E - mean_E);
				var_E /= n_replicas;
				double sigma_E = sqrt(var_E);
				T_min = (sigma_E > 1e-8) ? 1.0 / (2.0 * sigma_E) : 2.0;
				T_max = T_min * 10.0;
				temperatures.resize(n_replicas);
				for (int i = 0; i < n_replicas; ++i)
					temperatures[i] = T_min * pow(T_max / T_min, (double)i / (n_replicas - 1));
			}
			
			// Create replicas from bm actual
			std::vector<Rede> replicas(n_replicas, bm);  // Cada réplica é uma cópia da rede atual
			for (int i = 0; i < n_replicas; ++i)
				replicas[i].k = 1.0 / temperatures[i];
			
			auto start = high_resolution_clock::now();
			parallel_tempering_multi(
				replicas, bm_av_s, bm_av_ss,
				t_eq, t_step, relx, rept,
				n, mean, sigma, type, H,
				energy_per_replica, swap_acceptance_ratio, gen
			);
			auto end = high_resolution_clock::now();
			std::chrono::duration<double> duration_sec = end - start;
			times << duration_sec.count() << " " << swap_acceptance_ratio << std::endl;
	
			// --- Diagnóstico ---
			std::cout << "Swap acceptance rate: " << swap_acceptance_ratio << std::endl;
			for (int i = 0; i < n_replicas; ++i)
				std::cout << "T = " << temperatures[i] << ", <E> = " << energy_per_replica[i] << std::endl;
		}
	
		// Atualização dos parâmetros J e h
		for (int i = 0; i < bm.nbonds; ++i) {
			if (i < bm.n) {
				dh = eta_h * (bm_av_s[i] - av_s[i]);
				erroh += pow(bm_av_s[i] - av_s[i], 2);
				bm.h[i] -= dh;
			}
			dJ = eta_J * (bm_av_ss[i] - av_ss[i]);
			erroJ += pow(bm_av_ss[i] - av_ss[i], 2);
			bm.J[i] -= dJ;
		}
	
		erroJ = sqrt(erroJ / bm.nbonds);
		erroh = sqrt(erroh / bm.n);
		erros << inter << " " << setprecision(13) << erroJ << " " << setprecision(13) << erroh << endl;
	
		if (inter % cort == 0) {
			std::cout << text_name << " " << inter << " "
					  << "err_J" << " " << left << setw(13) << scientific << setprecision(6) << erroJ << " "
					  << "err_h" << " " << left << setw(13) << scientific << setprecision(6) << erroh << '\n';
		}
	
		if (erroJ <= min_erro_j && erroh <= min_erro_h) {
			cout << "Converged with: " << inter << " interations\n";
			break;
		}
	
		inter++;
	}
	
	// Close errors file
	erros.close();
	times.close();
	
//-----------------------------------------------------------------------------
	//Arquivo para salvar a rede obtida	
	ofstream network (file_rede_output.c_str());

	network << bm.n + 1 << endl;
	
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
	ofstream mag_corr (file_mag_corr_output.c_str());

	VecDoub bm_C(n*(n-1)/2, 0.0);

	//correlação pearson
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

	//Salva arquivo com Jij e Correlação-------------------
    //Nome do arquivo alvo
    

    //Abrindo arquivo output
    ofstream CorrJij (file_name_CorrJij.c_str());

	//Salvar arquivo com Jij e Pij
	//Nome do arquivo alvo
	
	//Abrindo arquivo output
	ofstream PJij (file_name_PJij.c_str());
	PJij << "Pij_ising" << endl;
    //Passando os valores para o arquivo
    for (int i = 0; i < bm.nbonds; i++)
    {
        CorrJij << left << setw(20) << setprecision(13) << bm.J[i] << 
                left << setw(18) << setprecision(13) << bm_C[i] << endl;

		PJij << bm.J[i] << " " << pearson_ising[i] << endl;
    }
	
//-----------------------------------------------------------------------------
//Salvar os dados separadamente

    //Arquivo para hi------------------
    

    //abrir arquivo
    ofstream file_hi (file_name_hi.c_str());

    //Arquivo para mi -----------------
    

    //abrir arquivo
    ofstream file_mi (file_name_mi.c_str());
	file_mi << "si_ising" << endl;
	file_hi << "hi_ising" << endl;
    for (int i = 0; i < n; i++)
    {
        file_hi << bm.h[i] << endl;
        file_mi << bm_av_s[i] << endl;
    }

    //fechando arquivos
    file_hi.close();
    file_mi.close();

    //-----------------------------------------------------

    //Arquivo para Jij------------------
    

    //abrir arquivo
    ofstream file_Jij (file_name_Jij.c_str());

    //Arquivo para Cij -----------------

    //abrir arquivo
    ofstream file_Cij (file_name_Cij.c_str());

	//Arquivo para Pij
	

	//abrir arquivo
	ofstream file_Pij (file_name_Pij.c_str());

	//Arquivo para sisj
	

	//abrir arquivo
	ofstream file_sisj (file_name_sisj.c_str());
    file_Jij << "Jij_ising" << endl;
	file_Cij << "Cij_ising" << endl;
	file_Pij << "Pij_ising" << endl;
	file_sisj << "sisj_ising" << endl;
	//salvar arquivos
    for (int i = 0; i < (n*(n-1)/2); i++)
    {
        file_Jij  << bm.J[i] << endl;
        file_Cij  << bm_C[i] << endl;
		file_Pij  << pearson_ising[i] << endl;
		file_sisj << bm_av_ss[i] << endl;
    }
	
    //fechando arquivos
    file_Jij.close();
    file_Cij.close();
	file_Pij.close();
    file_sisj.close();

	//-----------------------------------------------------
	//Calcular os tripletos

	int n_triplet = n*(n-1)*(n-2)/6;
	vector<double> bm_av_sss(n_triplet), ising_Triplet(n_triplet);
	vector<vector<double>> ising_M_av_ss(n, vector<double>(n));

	//Arquivo para Tijk
	

	//abrir arquivo
	ofstream file_Tijk (file_name_Tijk.c_str());

	//Arquivo para sisjsk
	

	//abrir arquivo
	ofstream file_sisjsk (file_name_sisjsk.c_str());

	file_Tijk << "Tijk_ising" << endl;
	file_sisjsk << "sisjsk_ising" << endl;
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
				
				file_Tijk   << ising_Triplet[ind] << endl;
				file_sisjsk << bm_av_sss[ind] << endl;
				
				ind++;
			}
		}
	}

	
	file_Tijk.close();
    file_sisjsk.close();

	return 0;
}


