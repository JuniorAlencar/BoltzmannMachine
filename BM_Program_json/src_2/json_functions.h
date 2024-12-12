#ifndef WRITE_JSON_H
#define WRITE_JSON_H


#include "network.h"
#include "create_folders.h"
#include <filesystem>
#include <nlohmann/json.hpp> // Inclua a biblioteca JSON antes das outras

using json = nlohmann::json;

struct exp_means {
    std::vector<double> av_s;
    std::vector<double> av_ss;
    std::vector<double> av_sss;
    std::vector<double> Cij_exp;
    std::vector<double> Pij_exp;
    std::vector<double> Tijk_exp;
};

class js_funct{
    public:
        inline void create_json_exp(const exp_means &data, const string &file_means);
        inline void write_json_properties(const string &filename,                 // Sample name
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
                            );
        exp_means load_json_exp(const std::string &file_means);
};


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
void js_funct::create_json_exp(const exp_means &data, const string &file_means){
    json j;
    j["S_exp"] = data.av_s;
    j["SS_exp"] = data.av_ss;
    j["SSS_exp"] = data.av_sss;
    j["Pij_exp"] = data.Pij_exp;
    j["Cij_exp"] = data.Cij_exp;
    j["Tijk_exp"] = data.Tijk_exp;
    
    std::ofstream file;
    file.open(file_means);
    file << j.dump(4);
    file.close();
};

inline exp_means js_funct::load_json_exp(const std::string &file_means){    
    // Abrindo o arquivo
    std::ifstream file(file_means);
    if (!file.is_open()) {
        throw std::runtime_error("Não foi possível abrir o arquivo: " + file_means);
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



class exp_mean_calculate{
    private:
        string use_exact;             // method exact/metropolis
        string filename;
    
    public:
        inline exp_means exp_calculate(const string &sample_input);
};

exp_means exp_mean_calculate::exp_calculate(const string &sample_input) {
    cout << "running exp_means" << endl;
	// Add flag class to use all functions defined
    js_funct rnd;
    
    // Using std::filesystem for path manipulation
    std::filesystem::path filePath(sample_input);
    
    // Remove the file extension
    std::filesystem::path withoutExtension = filePath.stem(); // Gets the filename without extension
    
    // Get the parent directory and combine with the file name without the extension
    std::filesystem::path resultPath = filePath.parent_path() / withoutExtension;

    // Convert the resulting path to a string
    std::string result = resultPath.string();
    
    // Check if file exp exists
    string file_means = result + ".json";

    // load struct with means exp values
    exp_means my_means;

    // Imprimir o caminho antes de verificar sua existência
    std::cout << "Verificando existência do arquivo: " << file_means << std::endl;

    fs::path file_path;
    
    if (!file_means.empty()) {
        file_path = fs::absolute(file_means);  // Ensure the path is valid and absolute
    } else {
        throw std::runtime_error("O caminho do arquivo está vazio.");
    }

    // Remover verificação de root path, apenas verificar se o caminho não está vazio
    if (file_means.empty()) {
        std::cerr << "Erro: Caminho inválido ou corrompido." << std::endl;
        return my_means;  // Retorna uma estrutura vazia
    }
    
    // Capturar qualquer exceção ao verificar a existência do arquivo
    try {
        if (fs::exists(file_path)) {
            cout << "Arquivo existe, abrindo..." << endl;
            // Load the experimental data from JSON file
            my_means = rnd.load_json_exp(file_means);
        } else {
            cout << "Arquivo experimental não existe, calculando..." << endl;
            
            string file_input = filePath.stem().string();
            
            ifstream data_input(sample_input.c_str());

            if (!data_input.is_open()) {
                cerr << "Erro ao abrir o arquivo: " << file_input << endl;
                return my_means;
            }

            string first_line, a;
            getline(data_input, first_line);
            int N = count(first_line.begin(), first_line.end(), ',');

            int nspins = N + 1;
            int m = 0;

            while (getline(data_input, first_line)) {
                m++;
            }
            m++;

            cout << "Número de spins = " << nspins << endl;
            cout << "Número de amostras = " << m << endl;

            data_input.clear();
            data_input.seekg(0, ios::beg);

            vector<vector<int>> M(m, vector<int>(N));

            for (int w = 0; w < m; w++) {
                for (int p = 0; p < N; p++) {
                    if (p < N - 1) {
                        getline(data_input, a, ',');
                    } else {
                        getline(data_input, a);
                    }
                    M[w][p] = stoi(a);
                }
            }

            data_input.close();

            int n_triplet = N * (N - 1) * (N - 2) / 6;
            int n_duplet = N * (N - 1) / 2;

            vector<double> s(N, 0.0), ss(n_duplet, 0.0), sss(n_triplet, 0.0);
            vector<vector<double>> M_ss(N, vector<double>(N));

            for (int p = 0; p < N; p++) {
                for (int w = 0; w < m; w++) {
                    s[p] += M[w][p];
                }
                s[p] /= m;
            }

            int ind = 0;
            for (int p = 0; p < N - 1; p++) {
                for (int pp = p + 1; pp < N; pp++) {
                    for (int w = 0; w < m; w++) {
                        ss[ind] += M[w][p] * M[w][pp];
                    }
                    ss[ind] /= m;
                    M_ss[p][pp] = ss[ind];
                    M_ss[pp][p] = M_ss[p][pp];
                    ind++;
                }
            }

            ind = 0;
            for (int i = 0; i < N - 2; i++) {
                for (int j = i + 1; j < N - 1; j++) {
                    for (int k = j + 1; k < N; k++) {
                        for (int s = 0; s < m; s++) {
                            sss[ind] += M[s][i] * M[s][j] * M[s][k];
                        }
                        sss[ind] /= m;
                        ind++;
                    }
                }
            }

            vector<double> Triplet(n_triplet, 0.0);
            ind = 0;
            for (int i = 0; i < N - 2; i++) {
                for (int j = i + 1; j < N - 1; j++) {
                    for (int k = j + 1; k < N; k++) {
                        Triplet[ind] = sss[ind] - s[i] * M_ss[j][k] - s[j] * M_ss[i][k]
                                       - s[k] * M_ss[i][j] + 2 * s[i] * s[j] * s[k];
                        ind++;
                    }
                }
            }

            vector<double> C(N * (N - 1) / 2, 0.0);
            vector<double> Pearson(n_duplet);
            ind = 0;
            for (int p = 0; p < N - 1; p++) {
                for (int pp = p + 1; pp < N; pp++) {
                    C[ind] = ss[ind] - s[p] * s[pp];
                    Pearson[ind] = C[ind] / sqrt((1 - pow(s[p], 2)) * (1 - pow(s[pp], 2)));
                    ind++;
                }
            }

            my_means.av_s = s;
            my_means.av_ss = ss;
            my_means.av_sss = sss;
            my_means.Cij_exp = C;
            my_means.Pij_exp = Pearson;
            my_means.Tijk_exp = Triplet;
            
            // Save experimental means in .json
            rnd.create_json_exp(my_means, file_means);
            
            // save mag_corr in .dat
            string file_name_output = result + "_mag_corr.dat";
            ofstream mag_corr(file_name_output.c_str());
            mag_corr << N << endl;

            for (int i = 0; i < N * (N - 1) / 2; i++) {
                mag_corr << " " << ss[i] << " " << C[i];
                if (i < N) {
                    mag_corr << " " << s[i];
                }
                mag_corr << endl;
            }

            mag_corr.close();
        }
    } catch (const fs::filesystem_error& ex) {
        std::cerr << "Erro ao verificar ou criar o arquivo: " << ex.what() << std::endl;
        return my_means;
    }

    return my_means;
};

#endif // !WRITE_JSON_H
