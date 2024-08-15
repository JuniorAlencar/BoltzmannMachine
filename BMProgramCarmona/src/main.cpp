#include <ctime>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <vector>
#include <boost/filesystem.hpp>

#include "exact_solution.h"
#include "experimental_means.h"
#include "forwardmethod.h"
#include "logger.h"
#include "metropolis.h"
#include "network.h"
#include "nr3.h"
#include "read_input_json.h"
#include "read_json.h"
#include "write_json.h"
#include "create_folders.h"

using namespace std;
namespace fs = boost::filesystem;

int main(int argc, char *argv[]) {
  // Create Results Folders
  create_folders();
  
  spdlog::flush_every(std::chrono::seconds(10));

  auto logger = get_logger();

  srand(time(NULL));
  string input_file = argv[1];
  // Read input variables in json file
  Params p;
  read_input_json(input_file, p);
  p.log_info();
  
  double beta = 1.0;
  
  // Matrix with all data (nspins X nrecords)
  vector<vector<int>> M;
  
  // Number spins
  int nspins;
  // Number samples
  int nrecords;
  
  read_tidy_data(p.input_data, M, nrecords, nspins);
  
  // Number of combinations in tuples (i,j) of system
  int nbonds = nspins * (nspins - 1) / 2;
  // Number of combinations in triples (i, j, k) of system
  int ntriplets = nspins * (nspins - 1) * (nspins - 2) / 6;
  
  // Carry first moment exp for each spin -> sigma_i_exp
  vector<double> av_s;
  // Carry second moment exp for each spin -> sigma_i*sigma_j_exp
  vector<double> av_ss;
  // Carry third moment exp for each spin _ sigma_i*sigma_j*sigma_k_exp
  vector<double> av_sss;

  // Carry first moment ising for each spin -> sigma_i_ising
  VecDoub bm_av_s(nspins, 0.0);
  // Carry second moment ising for each spin -> sigma_i*sigma_j_ising
  VecDoub bm_av_ss(nbonds, 0.0);
  // Carry third moment ising for each spin -> sigma_i*sigma_j*sigma_k_ising
  VecDoub_IO bm_av_sss(ntriplets);
  
  // Compute first, second, third moments, Cij, Pij, Triplet (Tijk) exp
  compute_exp_means(M, nrecords, nspins, av_s, av_ss, av_sss);
  
  // Network starts from zero with nspins
  Rede bm(nspins, 0, 0, 0, 0, 0);

  if (p.input_init_guess.find(".json") != string::npos) {
    bm = read_json(p.input_init_guess);
  }

  // MC variables
  int t_eq = p.t_eq;
  int t_rep = 2 * nspins;
  int n_rep = p.n_rep;
  int t_meas = p.t_meas;

  double err_SS = 1, err_S = 1;
  double dJ, dh;
  
  // Step to show results in display
  int iter_display = p.iter_display;
  
  // Interation initial
  int iter = 1;
  
  // Max number of interations
  int iter_max = p.iter_max;

  double eta_J, eta_J0 = p.eta_J;
  double eta_h, eta_h0 = p.eta_h;
  // min erro J
  double min_err_SS = p.min_err_SS;
  // min erro h
  double min_err_S = p.min_err_S;

  // Arquivo para salvar os erros ao longo do tempo

  ofstream cerr_f(p.output_err.c_str());

  // erros.seekg(0, std::ios_base::end);

  double av_av_s = 0.0;
  
  for (int i = 0; i < bm.n; ++i)
    av_av_s += av_s[i];
  
  av_av_s = abs(av_av_s / bm.n);

  double av_av_ss = 0.0;
  for (int i = 0; i < bm.nbonds; ++i)
    av_av_ss += av_ss[i];
  
  av_av_s = abs(av_av_ss / bm.nbonds);

  if (p.use_exact)
    logger->info("main: using full ensemble");
  else
    logger->info("main: using metropolis");
  
  while ((err_SS > min_err_SS || err_S > min_err_S) &&
         iter <= iter_max) //(inter <= inter_max)
  {
    srand(time(NULL) * time(NULL));

    //eta_J = eta_J0 * pow(iter, -0.4);
    //eta_h = eta_h0 * pow(iter, -0.4);
    err_SS = err_S = 0;
    eta_J = pow(iter, -0.4);
    eta_h = 2 * pow(iter, -0.4);

    if (p.use_exact) {
      exact_solution_bm(bm, bm_av_s, bm_av_ss, beta);
    } else {
      metropolis(bm, bm_av_s, bm_av_ss, bm_av_ss, t_eq, t_meas, t_rep, n_rep,
                 beta, false);
    }

    for (int i = 0; i < bm.nbonds; i++) {
      if (i < bm.n) {
        dh = eta_h * (bm_av_s[i] - av_s[i]);
        err_S += pow(bm_av_s[i] - av_s[i], 2);
        bm.h[i] -= dh;
      }

      dJ = eta_J * (bm_av_ss[i] - av_ss[i]);
      err_SS += pow(bm_av_ss[i] - av_ss[i], 2);
      bm.J[i] -= dJ;
    }

    err_SS = sqrt(err_SS / bm.nbonds);
    err_S = sqrt(err_S / bm.n);

    // Save erros
    cerr_f << std::setprecision(13) << iter << " " << err_SS  << " " << err_S << endl;

    if (iter % iter_display == 0 ||
        (err_SS < min_err_SS && err_S < min_err_S)) {
          
      logger->info("main: {0} iter={1:d} err_S={2:.6E} err_SS={3:.6E}",
                   p.run_name, iter, err_S, err_SS);
    }

    iter++;
  }

  // Fechar arquivos dos erros salvos
  cerr_f.close();


  if (p.use_exact) {
    exact_solution_triplet(bm, bm_av_s, bm_av_ss, bm_av_sss, beta);
  } else {
    metropolis(bm, bm_av_s, bm_av_ss, bm_av_sss, t_eq, t_meas, t_rep, n_rep,
               beta, true);
  }

  write_json(p.output_data, bm, av_s, av_ss, av_sss, bm_av_s, bm_av_ss,
             bm_av_sss);
  
  return 0;
}
