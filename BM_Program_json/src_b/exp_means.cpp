#include "exp_means.h"

exp_means exp_mean_calculate::exp_calculate(const string &filename) {
    cout << "running exp_means" << endl;
	// Add flag class to use all functions defined
    js_funct rnd;

    // Check if file exp exists
    string file_means = "../Results/" + filename + "_exp.json";

    // load struct with means exp values
    exp_means my_means;

    // Imprimir o caminho antes de verificar sua existência
    std::cout << "Verificando existência do arquivo: " << file_means << std::endl;

    fs::path file_path(file_means);

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

            string file_input = "../Data/TidyData/" + filename + ".dat";
            ifstream data_input(file_input.c_str());

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

            rnd.create_json_exp(my_means, filename);

            string file_name_output = "../Results/" + filename + "_mag_corr.dat";
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
