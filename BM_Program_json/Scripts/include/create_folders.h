#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <fmt/core.h>
#include <boost/filesystem.hpp>

using namespace std;
namespace fs = boost::filesystem;

// variable type is an auxiliary variable, where
// type=0 returns nothing,
// type=1 returns properties folder
// type=2 returns errors folder
// type=3 returns network folder
// type=4 returns specific heat folder
std::string create_folders(string &text_name, int &multiply_teq, int &multiply_relx, const bool &method, const int &type) {
	    
    // Count number of spins in sample ---------------------------------
    string file_input = "../Data/TidyData/" + text_name + ".dat";
    ifstream data_input (file_input.c_str());
    string first_line;
    string a;

    getline(data_input, first_line);
    int N = count(first_line.begin(), first_line.end(), ',');
    
    // Number of spins
    int nspins = N + 1; 
    
    // Converter para string ------------------------------------------	
	ostringstream os_teq;
	os_teq << multiply_teq*nspins;

	ostringstream os_relx;
	os_relx << multiply_relx*nspins;
	
	string teq_str = os_teq.str();
	string relx_str = os_relx.str();
    
    
    // Folders exact_solutions-----------------------------------------
    // Results
    // |
    // -----teq_{teq_value}
    //      |
    //       -----relx_{relx_value}
    //           |
    //            -----Properties, Specific_heat, network and errors
    // Set folder names using std::string
    string results_folder = "../Results";
    string teq_folder = results_folder + "/teq_" + teq_str;
    string relx_folder = teq_folder + "/relx_" + relx_str;
    string prop_folder = relx_folder + "/properties";
    string specific_heat_folder = relx_folder + "/specific_heat";
    string network_folder = relx_folder + "/network";
    string errors_folder = relx_folder + "/errors";

    // Create folders
    fs::create_directories(results_folder);
    fs::create_directories(teq_folder);
    fs::create_directories(relx_folder);
    fs::create_directories(prop_folder);
    fs::create_directories(specific_heat_folder);
    fs::create_directories(network_folder);
    fs::create_directories(errors_folder);

    // Create Data folders
    string Mag_Corr_folder = "../Data/Mag_Corr";
    fs::create_directories(Mag_Corr_folder);

    // Folders metropolis_solutions-----------------------------------------

    // Set folder names using std::string
    string results_folder_metropolis = "../Results_Metropolis";
    string teq_folder_metropolis = results_folder_metropolis + "/teq_" + teq_str;
    string relx_folder_metropolis = teq_folder_metropolis + "/relx_" + relx_str;
    string prop_folder_metropolis = relx_folder_metropolis + "/properties";
    string specific_heat_folder_metropolis = relx_folder_metropolis + "/specific_heat";
    string network_folder_metropolis = relx_folder_metropolis + "/network";
    string errors_folder_metropolis = relx_folder_metropolis + "/errors";

    // Create folders
    fs::create_directories(results_folder_metropolis);
    fs::create_directories(teq_folder_metropolis);
    fs::create_directories(relx_folder_metropolis);
    fs::create_directories(prop_folder_metropolis);
    fs::create_directories(specific_heat_folder_metropolis);
    fs::create_directories(errors_folder_metropolis);
    fs::create_directories(network_folder_metropolis);
    
    if(type == 0)
        return "";
    else if(type == 1)
        return method ? prop_folder : prop_folder_metropolis;
    else if(type == 2)
        return method ? errors_folder : errors_folder_metropolis;
    else if(type == 3)
        return method ? network_folder : network_folder_metropolis;
    else if(type == 4)
        return method ? specific_heat_folder : specific_heat_folder_metropolis;
    else
        return "enter with type value accept (int type <= 4)"

}