#ifndef READ_JSON_HPP
#define READ_JSON_HPP

#include "exp_means.h"  // Inclua o header com a definição da struct
#include "nr3.h"
#include <iostream>
#include <nlohmann/json.hpp>
#include "create_folders.hpp"
#include "network.h"
#include <stdexcept>
#include <fstream>

using json = nlohmann::json;

// Função inline para carregar o JSON e preencher a struct, retornando ambos
inline std::pair<exp_means, json> load_json_exp(const std::string &filename) {
    // Caminho completo do arquivo
    std::string full_filename = "../Results/" + filename + "_exp.json";

    // Abrindo o arquivo
    std::ifstream file(full_filename);
    if (!file.is_open()) {
        throw std::runtime_error("Não foi possível abrir o arquivo: " + full_filename);
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
    return std::make_pair(data, j);
}

#endif // !READ_JSON_HPP