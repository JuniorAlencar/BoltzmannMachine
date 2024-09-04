#ifndef WRITE_JSON_H
#define WRITE_JSON_H


#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <string>

#include "exp_means.h"
#include "network.h"
#include "nr3.h"

using namespace std;
using json = nlohmann::json;

// void write_json_properties(const std::string &filename,
//                 const int &multiply_relx,           // relx = multiply_relx*nspins
//                 const int &multiply_teq,            // t_eq=multiply_teq*nspins
//                 const int &nspins,                  // Number spins
//                 const bool &use_exact,              // If false==>Metropolis
//                 const double &min_err_h,            // Minimum value h to converge
//                 const double &min_err_j,            // Minimum value J to converge
//                 const Rede &bm,                     // Network
//                 const vector<double> &av_s,         // First moment exp
//                 const vector<double> &av_ss,        // Second moment exp
//                 const vector<double> &av_sss,       // Third moment exp
//                 const vector<double> &bm_av_s,      // First moment ising
//                 const vector<double> &bm_av_ss,     // Second moment ising
//                 const vector<double> &bm_av_sss,    // Third moment ising
//                 const vector<double> &C_exp,        // Covariance_exp
//                 const vector<double> &Pij_exp,      // Correlation_exp
//                 const vector<double> &Tijk_exp,     // Triplet_exp
//                 const vector<double> &C_ising,      // Covariance_ising
//                 const vector<double> &Pij_ising,    // Correlation_ising
//                 const vector<double> &Tijk_ising)  // Triplet_ising 

// filename: sample name
// multiply_relx: relx=multiply_relx*nspins
// multiply_t_eq: t_eq=multiply_t_eq*nspins
// nspins: number of spins in sample
// use_exact: true -> exact_solutions, false -> metropolis
// min_err_h: minimum error in h to MC converge
// min_err_J: minimum error in J to MC converge
// exp_means: allocate all means experimentals (see exp_means.h)
// bm_av_s: First moment Ising
// bm_av_ss: Second moment Ising
// bm_av_sss: Third moment Ising
// C_ising: Covariance ising
// Pij_ising: Correlation ising (Pearson)
// Tijk: Triplet Ising
void write_json_properties(const std::string &filename,
                const int &relx,           // relx = multiply_relx*nspins
                const int &teq,            // t_eq=multiply_teq*nspins
                const int &nspins,                  // Number spins
                const bool &use_exact,              // If false==>Metropolis
                const double &min_err_h,            // Minimum value h to converge
                const double &min_err_j,            // Minimum value J to converge
                const Rede &bm,
                const exp_means &Means_Exp,         // Means_Experimental           
                const VecDoub_IO &bm_av_s,      // First moment ising
                const VecDoub_IO &bm_av_ss,     // Second moment ising
                const vector<double> &bm_av_sss,    // Third moment ising
                const VecDoub &C_ising,      // Covariance_ising
                const vector<double> &Pij_ising,    // Correlation_ising
                const vector<double> &Tijk_ising   // Triplet_ising 
                )  
{
    int n=bm.n;
    int npairs=bm.nbonds;
    int ntriplets=bm_av_sss.size();
    
    vector<int> s(n,0);
    vector<double> bmavs(n,0), bmavss(npairs,0);
    vector<double> h(n,0.0), J(npairs,0.0);
    vector<double> C_is(npairs, 0.0);

    // h_i and First Moment
    for (int i = 0; i < n; ++i){
        s[i] = bm.s[i];
        bmavs[i] = bm_av_s[i];
    }
    
    // Second moment and J_ij
    for (int i = 0; i < npairs; ++i){
        bmavss[i] = bm_av_ss[i];
        C_is[i] = C_ising[i];
    }
    
    // Converter to string mantendo notação cientifica
	ostringstream os_j;
    os_j << scientific << setprecision(2) << min_err_j;

	ostringstream os_h;
    os_h << scientific << setprecision(2) << min_err_h;

    string min_erro_j_str = os_j.str();
	string min_erro_h_str = os_h.str();

    json j;
 
    // Number of spins in data
    j["n_spins"] = nspins;
    
    // First moment exp
    j["S_exp"] = Means_Exp.bm_av_s;
    // Second moment exp
    j["SS_exp"] = Means_Exp.bm_av_ss;
    // Third moment exp
    j["SSS_exp"] = Means_Exp.bm_av_sss;
    
    // First moment Ising
    j["S_ising"] = bmavs;
    // Second moment Ising
    j["SS_ising"] = bmavss;
    // Third moment Ising
    j["SSS_ising"] = bm_av_sss;
    
    // Covariance exp
    j["Cov_exp"] = Means_Exp.C_exp;
    // Correlation exp
    j["Corr_exp"] = Means_Exp.Pij_exp;
    // Triplet exp
    j["Tijk_exp"] = Means_Exp.Tijk_exp;
    
    // Covariance ising
    j["Cov_ising"] = C_is;
    // Correlation ising
    j["Corr_ising"] = Pij_ising;
    // Triplet ising
    j["Tijk_ising"] = Tijk_ising;
    
    string file_name;
    
    string teq_str = to_string(teq);
    string relx_str = to_string(relx);
    
    if(use_exact == true)
        file_name = "../Results/" + filename + "/teq_" + teq_str + 
                    "/relx_" + relx_str + "/system_analysis/properties/props_min_h_" +
                    min_erro_h_str + "_min_j_" + min_erro_j_str + ".json";
    else
        file_name = "../Results_Metropolis/" + filename + "/teq_" + teq_str + 
            "/relx_" + relx_str + "/system_analysis/properties/props_min_h_" +
            min_erro_h_str + "_min_j_" + min_erro_j_str + ".json";
    
    ofstream file;
    
    file.open(file_name);
    file << j.dump(4);  // The '4' here is for pretty-printing with an indent of 4 spaces
    file.close();

}

#endif // WRITE_JSON_H