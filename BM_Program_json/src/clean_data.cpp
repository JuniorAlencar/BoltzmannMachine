#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>

using namespace std;

int main(int argc, char *argv[]){
    string file_input	= argv[1];
    
    size_t pos_first = file_input.find_first_of("/\\");
    size_t pos_second = file_input.find_first_of("/\\", pos_first + 1);
    size_t pos_third = file_input.find_first_of("/\\", pos_second + 1);
    // Obter o caminho a partir da segunda barra
    string path_from_third = file_input.substr(pos_third + 1);
    
    string file_output = "../Results/" + path_from_third;
    
    ifstream data_input(file_input.c_str());
    
    if (!data_input.is_open()) {
        cerr << "Erro ao abrir o arquivo de entrada!" << endl;
        return 1;
    }
    
    ofstream data_output(file_output.c_str());
    
    if (!data_output.is_open()) {
        cerr << "Erro ao criar o arquivo de saída!" << endl;
        return 1;
    }

    string line;
    string a;

    // Ler o arquivo linha por linha e processar
    while (getline(data_input, line))
    {
        stringstream ss(line);  // Transformar a linha em stream para processar os valores
        bool first = true;  // Controlar o primeiro valor para não adicionar espaço antes dele
        
        while (getline(ss, a, ' '))  // Usar espaço como delimitador
        {
            double value = stod(a);  // Converter string para double
            
            // Condicional: Se o valor for menor que 0.5, coloca +1; caso contrário, -1
            if (!first) {
                data_output << ",";  // Adicionar espaço antes de cada valor, exceto o primeiro
            }
            first = false;

            if (value < 0.5) {
                data_output << "1";
            } else {
                data_output << "-1";
            }
        }
        data_output << "\n";  // Nova linha após cada linha da matriz
    }

    // Fechar os arquivos
    data_input.close();
    data_output.close();

    cout << "Matriz salva em " << file_output << endl;

    return 0;
}
