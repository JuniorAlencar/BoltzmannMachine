#ifndef EXP_MEANS_H
#define EXP_MEANS_H

#include <string>
#include <vector>

void read_tidy_data(const std::string &filename, 
                    std::vector<std::vector<int>> &M,
                    int &nrows, int &ncols);
void compute_exp_means(const std::vector<std::vector<int>> &M,
                       const int &nrows,
                       const int &ncols,
                       std::vector<double> &s, 
                       std::vector<double> &ss,
                       std::vector<double> &sss);
#endif