#include "json_functions.h"

void js_funct::write_json_properties( const string &filename,                 // Sample name
                            const int &nspins,
                            const int &n_duplet,                         // Network
                            // MC parameters
                            const int &multi_relx,                        // Parameter to MC
                            const int &multi_teq,                         // Parameter to MC
                            const double &jmin,                     // Parameter to MC
                            const double &hmin,                     // Parameter to MC
                            const string &method,                   // exact or metropolis
                            // Exp properties
                            const exp_means &data_exp,              //experimental means
                            // Ising properties
                            const VecDoub &bm_av_s,                 // First moment Ising
                            const VecDoub &bm_av_ss,                // Second moment Ising
                            const std::vector<double> &bm_av_sss,   // Third moment Ising
                            const std::vector<double> &Pij_ising,   // Correlation Ising
                            const std::vector<double> &Tijk_ising   // Triplet Ising
                            ){
    // ==> The json library accepts vector<double> in arguments <==
    // Converter bm_av_s and bm_av_ss from VecDoub to vector<double>
    vector<double> bmavs(nspins,0), bmavss(n_duplet,0);
    // Converter Cij_ising from VecDoub to vector<double>
    vector<double> C_ij_ising(n_duplet,0);
    
    // Converting VecDoub type vectors to vector<double>-----------------
    for(int i = 0; i < nspins; ++i)
        bmavs[i] = bm_av_s[i];
    
    for(int i = 0; i<n_duplet; ++i)
        bmavss[i] = bm_av_ss[i];
    
    int aux = 0;
    for(int i = 0; i < nspins-1; ++i ){
        for(int j = i+1; j < nspins; ++j){
            C_ij_ising[aux] = bmavss[aux] - bmavs[i]*bmavs[j];
            aux ++;
        }
    }
    
    json j;
    j["sample"] = filename;
    j["n_spins"] = nspins;
    // Monte Carlo parameters -----------------
    j["MC"] = json::object();  // Initialize "MC" as an empty object
    j["MC"]["relx"] = multi_relx*nspins;
    j["MC"]["teq"] = multi_teq*nspins;
    j["MC"]["jmin"] = jmin;
    j["MC"]["hmin"] = hmin;
    
    // Method used -----------------
    j["method"] = method;
    // First, second and third moment

    j["S_exp"] = data_exp.av_s;

    j["SS_exp"] = data_exp.av_ss;

    j["SSS_exp"] = data_exp.av_sss;

    // Correlation, covariance and triplet

    j["Pij_exp"] = data_exp.Pij_exp;

    j["Cij_exp"] = data_exp.Cij_exp;

    j["Tijk_exp"] = data_exp.Tijk_exp;

    // Ising means -----------------------------
    // First, second and third moment
    j["S_ising"] = bmavs;
    j["SS_ising"] = bmavss;
    j["SSS_ising"] = bm_av_sss;
    // Correlation, covariance and triplet
    j["Pij_ising"] = Pij_ising;
    j["Cij_ising"] = C_ij_ising;
    j["Tijk_ising"] = Tijk_ising;
    
    // Return repository to properties
    int type = 1;
    c_folders rnd;
    string full = rnd.create_folders(filename, multi_teq, multi_relx, method, type);
    string full_filename = full + "/" + filename + ".json";
    
    std::ofstream file;
    file.open(full_filename);
    file << j.dump(4);
    file.close();
};

// Function to save exp_means values
void js_funct::create_json_exp(const exp_means &data, const string &filename){
    json j;
    j["S_exp"] = data.av_s;
    j["SS_exp"] = data.av_ss;
    j["SSS_exp"] = data.av_sss;
    j["Pij_exp"] = data.Pij_exp;
    j["Cij_exp"] = data.Cij_exp;
    j["Tijk_exp"] = data.Tijk_exp;
    // save files exp
    string full_filename = "../Results/" + filename + "_exp.json"; 
    
    std::ofstream file;
    file.open(full_filename);
    file << j.dump(4);
    file.close();
};

exp_means js_funct::load_json_exp(const std::string &filename){    
    // Abrindo o arquivo
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Não foi possível abrir o arquivo: " + filename);
    }

    // Carregando o conteúdo do arquivo em um objeto JSON
    json j;
    file >> j;

    // Criando a struct exp_means e atribuindo os valores do JSON
    exp_means data;
    data.av_s = j["S_exp"].get<std::vector<double>>();
    data.av_ss = j["SS_exp"].get<std::vector<double>>();
    data.av_sss = j["SSS_exp"].get<std::vector<double>>();
    data.Pij_exp = j["Pij_exp"].get<std::vector<double>>();
    data.Cij_exp = j["Cij_exp"].get<std::vector<double>>();
    data.Tijk_exp = j["Tijk_exp"].get<std::vector<double>>();

    // Fechando o arquivo
    file.close();

    // Retornando a struct exp_means e o próprio JSON
    return data;
    };