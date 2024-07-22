#include <ctime>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

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

using namespace std;

int main(int argc, char *argv[]) {
  spdlog::flush_every(std::chrono::seconds(10));

  auto logger = get_logger();

  srand(time(NULL));
  string input_file = argv[1];

  Params p;
  read_input_json(input_file, p);
  p.log_info();
  double beta = 1.0;
  vector<vector<int>> M;
  int nspins, nrecords;
  read_tidy_data(p.input_data, M, nrecords, nspins);
  int nbonds = nspins * (nspins - 1) / 2;
  int ntriplets = nspins * (nspins - 1) * (nspins - 2) / 6;

  vector<double> av_s, av_ss, av_sss;
  VecDoub bm_av_s(nspins, 0.0), bm_av_ss(nbonds, 0.0);
  VecDoub_IO bm_av_sss(ntriplets);

  compute_exp_means(M, nrecords, nspins, av_s, av_ss, av_sss);

  Rede bm(nspins, 0, 0, 0, 0, 0);

  if (p.input_init_guess.find(".json") != string::npos) {
    bm = read_json(p.input_init_guess);
  }

  // variaveis para MC
  int t_eq = p.t_eq;
  int t_rep = 2 * nspins;
  int n_rep = p.n_rep;
  int t_meas = p.t_meas;

  double err_SS = 1, err_S = 1;
  double dJ, dh;
  // int inter_ini = atoi(argv[3]);
  // int inter_max = atoi(argv[4]);
  int iter_display = p.iter_display;
  int iter = 1;              //= inter_ini;
  int iter_max = p.iter_max; // mesmo que meu step_relax

  double eta_J, eta_J0 = p.eta_J; // atof(argv[2]);
  double eta_h, eta_h0 = p.eta_h;

  double min_err_SS = p.min_err_SS;
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

    eta_J = eta_J0 * pow(iter, -0.4);
    eta_h = eta_h0 * pow(iter, -0.4);

    if (p.use_exact && nspins < 24) {
      exact_solution_bm(bm, bm_av_s, bm_av_ss, beta);
    } else {
      metropolis(bm, bm_av_s, bm_av_ss, bm_av_ss, t_eq, t_meas, t_rep, n_rep,
                 beta, false);
    }

    err_SS = err_S = 0;
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

    err_SS = sqrt(err_SS) / bm.nbonds / av_av_ss;
    err_S = sqrt(err_S) / bm.n / av_av_s;

    // Salvando Erros
    cerr_f << std::setprecision(13) << iter << "," << err_S << "," << err_SS
           << '\n';

    if (iter % iter_display == 0 ||
        (err_SS < min_err_SS && err_S < min_err_S)) {
      logger->info("main: {0} iter={1:d} err_S={2:07.6f} err_SS={3:07.6f}",
                   p.run_name, iter, err_S, err_SS);
    }

    iter++;
  }

  // Fechar arquivos dos erros salvos
  cerr_f.close();

  if (p.use_exact && nspins < 24) {
    exact_solution_triplet(bm, bm_av_s, bm_av_ss, bm_av_sss, beta);
  } else {
    metropolis(bm, bm_av_s, bm_av_ss, bm_av_sss, t_eq, t_meas, t_rep, n_rep,
               beta, true);
  }

  write_json(p.output_data, bm, av_s, av_ss, av_sss, bm_av_s, bm_av_ss,
             bm_av_sss);

  return 0;
}
