#include "create_folders.h"

std::string c_folders::create_folders(const string &text_name, const int &multiply_teq, const int &multiply_relx, const string &method, const int &type) {
    // Contar número de spins no arquivo de amostra
    string file_input = "../Data/TidyData/" + text_name + ".dat";
    ifstream data_input(file_input.c_str());

    if (!data_input.is_open()) {
        cerr << "Erro ao abrir o arquivo: " << file_input << endl;
        return "";
    }

    string first_line;
    getline(data_input, first_line);
    int N = count(first_line.begin(), first_line.end(), ',');
    int nspins = N + 1;

    // Converter valores de multiplicação para string
    ostringstream os_teq, os_relx;
    os_teq << multiply_teq * nspins;
    os_relx << multiply_relx * nspins;
    string teq_str = os_teq.str();
    string relx_str = os_relx.str();
    
    // Definir os nomes das pastas
    string Results_folder = "../Results";
    string results_folder = Results_folder + "/" + method;
    string teq_folder = results_folder + "/teq_" + teq_str;
    string relx_folder = teq_folder + "/relx_" + relx_str;
    string prop_folder = relx_folder + "/properties";
    string specific_heat_folder = relx_folder + "/specific_heat";
    string network_folder = relx_folder + "/network";
    string errors_folder = relx_folder + "/errors";

    // Se type for diferente de 0, não criar as pastas, apenas retornar o caminho correspondente
    if (type != 0) {
        switch (type) {
            case 1: return prop_folder;
            case 2: return errors_folder;
            case 3: return network_folder;
            case 4: return specific_heat_folder;
            case 5: return Results_folder;
            default: return "Entre com um valor de tipo válido (int type <= 5)";
        }
    }

    // Se type for 0, criar as pastas
    try {
        std::cout << "Verificando existência do diretório: " << Results_folder << std::endl;
        if (!fs::exists(Results_folder)) {
            std::cout << "Criando Results_folder: " << Results_folder << std::endl;
            fs::create_directory(Results_folder);  // Testar a criação do diretório
        }

        std::cout << "Criando results_folder: " << results_folder << std::endl;
        fs::create_directories(results_folder);  // Cria os demais diretórios

        std::cout << "Criando teq_folder: " << teq_folder << std::endl;
        fs::create_directories(teq_folder);

        std::cout << "Criando relx_folder: " << relx_folder << std::endl;
        fs::create_directories(relx_folder);

        std::cout << "Criando prop_folder: " << prop_folder << std::endl;
        fs::create_directories(prop_folder);

        std::cout << "Criando specific_heat_folder: " << specific_heat_folder << std::endl;
        fs::create_directories(specific_heat_folder);

        std::cout << "Criando network_folder: " << network_folder << std::endl;
        fs::create_directories(network_folder);

        std::cout << "Criando errors_folder: " << errors_folder << std::endl;
        fs::create_directories(errors_folder);
    } catch (const fs::filesystem_error& ex) {
        std::cerr << "Erro ao criar pastas em: " << ex.path1() << " -> " << ex.what() << std::endl;
        return "";
    }

    // Se type for 0, retornar vazio após a criação das pastas
    return "";
}
