#include "read_input_json.h"
#include <fstream>
#include <iostream>

using json = nlohmann::json;
using namespace std;

void read_input_json(const string &filename, Params &p)
{

    try
    {
        std::ifstream file(filename.c_str());
        json j;
        file >> j;
        file.close();

        p.run_name = j["run_name"];
        p.input_data = j["input_data"];
        p.input_init_guess = j["input_init_guess"];
        p.output_data = j["output_data"];
        p.output_err = j["output_err"];
        p.use_exact = j["use_exact"];

        p.iter_display = j["relax"]["iter_display"];
        p.iter_max = j["relax"]["iter_max"];
        p.eta_J = j["relax"]["eta_J"];
        p.eta_h = j["relax"]["eta_h"];
        p.min_err_S = j["relax"]["min_err_S"];
        p.min_err_SS = j["relax"]["min_err_SS"];
        p.n_rep = j["mc"]["n_rep"];
        p.t_eq = j["mc"]["t_eq"];
        p.t_meas = j["mc"]["t_meas"];

    }
    catch (const json::parse_error &e)
    {
        std::cerr << "JSON parsing error: " << e.what() << "\n";
        exit(1);
    }
    catch (const std::ifstream::failure &e)
    {
        std::cerr << "File I/O error: " << e.what() << "\n";
        exit(1);
    }
    catch (...)
    {
        std::cerr << "An unknown error occurred\n";
        exit(1);
    }
}