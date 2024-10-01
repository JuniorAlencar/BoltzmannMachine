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
std::string create_folders(const string &text_name,const int &multiply_teq, const int &multiply_relx, const bool &method, const int &type) {
	    
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
    // If method = exact, create results, else create results_metropolis
    string results_folder = method ? "../Results" : "../Results_Metropolis";
    
    string teq_folder = results_folder + "/teq_" + teq_str;
    string relx_folder = teq_folder + "/relx_" + relx_str;
    string prop_folder = relx_folder + "/properties";
    string specific_heat_folder = relx_folder + "/specific_heat";
    string network_folder = relx_folder + "/network";
    string errors_folder = relx_folder + "/errors";

    // Create Data folders
    string Mag_Corr_folder = "../Data/Mag_Corr";
    fs::create_directories(Mag_Corr_folder);
    // if type = 0, just create folders and return nothing
    if(type == 0){
         // Create folders
        fs::create_directories(results_folder);
        fs::create_directories(teq_folder);
        fs::create_directories(relx_folder);
        fs::create_directories(prop_folder);
        fs::create_directories(specific_heat_folder);
        fs::create_directories(network_folder);
        fs::create_directories(errors_folder);
        return "";}
    else if(type == 1)
        return prop_folder;
    else if(type == 2)
        return errors_folder;
    else if(type == 3)
        return network_folder;
    else if(type == 4)
        return specific_heat_folder;
    else
        return "enter with type value accept (int type <= 4)";
}