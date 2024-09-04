#ifndef CREATEFOLDERS_H
#define CREATEFOLDERS_H

#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;
namespace fs = boost::filesystem;

// Create folders with tree:
// Results/Results_Metropolis
// |
// ---teq_{nspins*multiply_teq}
//    |
//     ---relx_{nspins*multiply_relx}
//     |
//     ----Network
//     ----Erro
//     ----system_analysis (values properties in .json format)
void create_folders(const string &filename, const int &multiply_teq, 
                    const int &multiply_relx, const int &nspins, bool &use_exact) {
    
    int teq = multiply_teq * nspins;
    int relx = multiply_relx * nspins;

    // Create Data folders
    string Mag_Corr_folder;
    
    
    // Folders exact_solutions-----------------------------------------
    // Set folder names using std::string
    string results_folder;
 
    // Create folders to specific teq and relx
    string teq_str = to_string(teq);
    string relx_str = to_string(relx);
    
    
    if(use_exact == true){
        results_folder = "../Results";
        Mag_Corr_folder = "../Data/Mag_Corr/" + filename + "/exact/teq_" + teq_str + "/relx_" + relx_str;

    }
    
    else{
        results_folder = "../Results_Metropolis";
        Mag_Corr_folder = "../Data/Mag_Corr/" + filename + "/metropolis/teq_" + teq_str + "/relx_" + relx_str;
    }
    
    string results_file = results_folder + "/" + filename;
    string teq_folder = results_file + "/teq_" + teq_str;
    string relx_folder = teq_folder + "/relx_" + relx_str;
    string network_folder = relx_folder + "/Network";
    string erro_folder = relx_folder + "/Erro";
    string system_analysis_folder = relx_folder + "/system_analysis";
    string specific_heat_folder = system_analysis_folder + "/specific_heat";
    string properties_folder = system_analysis_folder + "/properties";

    // Create folders
    fs::create_directories(Mag_Corr_folder);
    fs::create_directories(results_folder);
    fs::create_directories(results_file);
    fs::create_directories(erro_folder);
    fs::create_directories(network_folder);
    fs::create_directories(erro_folder);
    fs::create_directories(specific_heat_folder);
    fs::create_directories(properties_folder);
}

#endif /* _CREATEFOLDERS_H_ */