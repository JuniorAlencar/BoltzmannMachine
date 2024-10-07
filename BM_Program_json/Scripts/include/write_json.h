#ifndef WRITE_JSON_H
#define WRITE_JSON_H

#include <vector>
#include <iostream>
#include <fstream> // Para manipulação de arquivos
#include <nlohmann/json.hpp> // Inclua a biblioteca JSON antes das outras
#include "exp_means.h"
#include "nr3.h"
#include "create_folders.h"
#include "network.h"

using json = nlohmann::json;

using namespace std;

// Function to save exp_means values
json create_json_exp(const exp_means &data, const string &filename) {
    json j;
    j["S_exp"] = data.av_s;
    j["SS_exp"] = data.av_ss;
    j["SSS_exp"] = data.av_sss;
    j["Pij_exp"] = data.Pij_exp;
    j["Cij_exp"] = data.Cij_exp;
    j["Tijk_exp"] = data.Tijk_exp;
    // save files exp
    string full_filename = "../Results/" + filename + "_exp.json"; 
    string full_filename_metro = "../Results_metropolis/" + filename + "_exp.json"; 
    
    std::ofstream file;
    file.open(full_filename);
    file << j.dump(4);
    file.close();

    std::ofstream file_metro;
    file_metro.open(full_filename_metro);
    file_metro << j.dump(4);
    file_metro.close();
    
    return j;
};

void write_json_properties( const string &filename,                 // Sample name
                            const Rede &bm,                         // Network
                            // MC parameters
                            const int &multi_relx,                        // Parameter to MC
                            const int &multi_teq,                         // Parameter to MC
                            const double &jmin,                     // Parameter to MC
                            const double &hmin,                     // Parameter to MC
                            const bool &method,                     // if True ->exact, else ->metropolis
                            // Exp properties
                            const std::vector<double> &av_s,        // First moment exp
                            const std::vector<double> &av_ss,       // Second moment exp
                            const std::vector<double> &av_sss,      // Third moment exp
                            const std::vector<double> &Cij_exp,     // Covariance exp
                            const std::vector<double> &Pij_exp,     // Correlation exp
                            const std::vector<double> &Tijk_exp,    // Triplet exp
                            // Ising properties
                            const VecDoub &bm_av_s,                 // First moment Ising
                            const VecDoub &bm_av_ss,                // Second moment Ising
                            const std::vector<double> &bm_av_sss,   // Third moment Ising
                            const VecDoub &Cij_ising,               // Covariance Ising
                            const std::vector<double> &Pij_ising,   // Correlation Ising
                            const std::vector<double> &Tijk_ising   // Triplet Ising
                            ){
    // ==> The json library accepts vector<double> in arguments <==
    int n=bm.n;
    int npairs=bm.nbonds;
    // Converter bm_av_s and bm_av_ss from VecDoub to vector<double>
    vector<double> bmavs(n,0), bmavss(npairs,0);
    // Converter Cij_ising from VecDoub to vector<double>
    vector<double> C_ij_ising(npairs,0);
    
    // Converting VecDoub type vectors to vector<double>-----------------
    for(int i = 0; i < n; ++i)
        bmavs[i] = bm_av_s[i];
    
    for(int i = 0; i<npairs; ++i)
        bmavss[i] = bm_av_ss[i];
    
    int aux = 0;
    for(int i = 0; i < n-1; ++i ){
        for(int j = i+1; j < n; ++j){
            C_ij_ising[aux] = bmavss[aux] - bmavs[i]*bmavs[j];
            aux ++;
        }
    }
    // If method = true, return 'exact', else return 'metropolis'
    string Method = method ? "exact" : "metropolis";
    
    json j;
    j["sample"] = filename;
    j["n_spins"] = n;
    // Monte Carlo parameters -----------------
    j["MC"] = json::object();  // Initialize "MC" as an empty object
    j["MC"]["relx"] = multi_relx*n;
    j["MC"]["teq"] = multi_teq*n;
    j["MC"]["jmin"] = jmin;
    j["MC"]["hmin"] = hmin;
    
    // Method used -----------------
    j["method"] = Method;
    // First, second and third moment

    j["S_exp"] = av_s;

    j["SS_exp"] = av_ss;

    j["SSS_exp"] = av_sss;

    // Correlation, covariance and triplet

    j["Pij_exp"] = Pij_exp;

    j["Cij_exp"] = Cij_exp;

    j["Tijk_exp"] = Tijk_exp;
    // Correlation, covariance and triplet
    j["Pij_exp"] = Pij_exp;
    j["Cij_exp"] = Cij_exp;
    j["Tijk_exp"] = Tijk_exp;

    // Ising means -----------------------------
    // First, second and third moment
    j["S_ising"] = bmavs;
    j["SS_ising"] = bmavss;
    j["SSS_ising"] = av_sss;
    // Correlation, covariance and triplet
    j["Pij_ising"] = Pij_ising;
    j["Cij_ising"] = C_ij_ising;
    j["Tijk_ising"] = Tijk_ising;
    
    // Return repository to properties
    int type = 1;
    
    string full = create_folders(filename, multi_teq, multi_relx, method, type);
    string full_filename = full + "/" + filename + ".json";
    
    std::ofstream file;
    file.open(full_filename);
    file << j.dump(4);
    file.close();
};

#endif // WRITE_JSON_H