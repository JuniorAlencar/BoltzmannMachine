#include "metropolis.h"

void metropolis(Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss,
                VecDoub_IO &av_sss, const int &t_eq, const int &t_meas,
                const int &t_rep, const int &n_rep, const double &beta,
                const bool &tr) {
  int s_flip;
  double p;
  double dE;

  // tem um micro estado inicial aleatório
  // int mmm = 0.0;
  // for (int j = 0; j < r.n; j++)
  // 	mmm += r.s[j];
  // std::cout << "mag total = " << mmm << std::endl;
  int t_step = t_eq + t_meas;

  for (int rep = 0; rep < n_rep; rep++) {
    for (int i = 0; i < r.n; ++i) {
      r.s[i] = -1;
      if ((double)rand() / RAND_MAX > 0.5)
        r.s[i] = 1;
    }

    // uma rodada de mc inteira - não gera um estado inicial novo...
    for (int t = 0; t < t_step; t++) {
      // i >= t_eq e de relx em relx acumula medias .... relx tem a ver
      // com o número de flips
      // acumulando primeiro e segundo momento
      if (t >= t_eq && (int)(t - t_eq) % (int)t_rep == 0) {
        for (int i = 0; i < r.n; i++)
          av_s[i] += r.s[i];

        int ind = 0;
        for (int i = 0; i < r.n - 1; i++) {
          for (int j = i + 1; j < r.n; j++) {
            av_ss[ind] += r.s[i] * r.s[j];
            ind++;
          }
        }

        if (tr == true) {
          ind = 0;
          for (int i = 0; i < r.n - 2; i++) {
            for (int j = i + 1; j < r.n - 1; j++) {
              for (int k = j + 1; k < r.n; k++) {
                av_sss[ind] += r.s[i] * r.s[j] * r.s[k];
                ind++;
              }
            }
          }
        }
      }

      s_flip = rand() % r.n; // flip aleatoriamente
      dE = delta_E(r, s_flip);

      if (dE <= 0) {
        r.s[s_flip] = -r.s[s_flip];
      } else {
        p = ((double)rand() / RAND_MAX);

        if (p < exp(-beta * dE)) {
          r.s[s_flip] = -r.s[s_flip];
        }
      }
    }
  }

  // t_step/relx = número de acumulações para cada repetição
  // rept*(t_step/relx) = número total de acumulações
  for (int i = 0; i < r.n; i++)
    av_s[i] /= (n_rep * t_meas / t_rep);

  int ind = 0;
  for (int i = 0; i < r.n - 1; i++) {
    for (int j = i + 1; j < r.n; j++) {
      av_ss[ind] /= (n_rep * t_meas / t_rep);
      ind++;
    }
  }

  if (tr == true) {
    ind = 0;
    for (int l = 0; l < r.n - 2; l++) {
      for (int j = l + 1; j < r.n - 1; j++) {
        for (int k = j + 1; k < r.n; k++) {
          av_sss[ind] /= (n_rep * t_meas / t_rep);
          ind++;
        }
      }
    }
  }
}
