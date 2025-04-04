#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;
namespace fs = boost::filesystem;

void create_folders(string method) {
    // Folders exact_solutions-----------------------------------------
    // Set folder names using std::string
    string tests_folder = "../tests"; 
    string results_folder = "../Results_" + method;
    string specificHeat_folder = results_folder + "/SpecificHeat";
    string comparative_folder = results_folder + "/Comparative";
    string CorrJij_folder = results_folder + "/CorrJij";
    string Energy_folder = results_folder + "/Energy";
    string Erro_folder = results_folder + "/Erro";
    string Mag_Corr_ising_folder = results_folder + "/Mag_Corr_ising";
    string Magnetization_vs_T_folder = results_folder + "/Magnetization_vs_T";
    string MatrixJij_folder = results_folder + "/MatrixJij";
    string Network_folder = results_folder + "/Network";
    string PJij_folder = results_folder + "/PJij";
    string SeparateData_folder = results_folder + "/SeparateData";

    // Create folders
    fs::create_directories(tests_folder);
    fs::create_directories(results_folder);
    fs::create_directories(specificHeat_folder);
    fs::create_directories(comparative_folder);
    fs::create_directories(CorrJij_folder);
    fs::create_directories(Energy_folder);
    fs::create_directories(Erro_folder);
    fs::create_directories(Mag_Corr_ising_folder);
    fs::create_directories(Magnetization_vs_T_folder);
    fs::create_directories(MatrixJij_folder);
    fs::create_directories(Network_folder);
    fs::create_directories(PJij_folder);
    fs::create_directories(SeparateData_folder);

    // folders method to tests
    string tests_method = tests_folder + "/" + method;
    fs::create_directories(tests_method);
    
    // Create subfolders
    string hi_tests = tests_method + "/hi";
    string Jij_tests = tests_method + "/Jij";
    string H_tests = tests_method + "/H";
    string si_tests = tests_method + "/si";
    string sisj_tests = tests_method + "/sisj";
    
    fs::create_directories(hi_tests);
    fs::create_directories(Jij_tests);
    fs::create_directories(si_tests);
    fs::create_directories(H_tests);
    fs::create_directories(sisj_tests);

    // Add more folders as required
    string correlation_folder = comparative_folder + "/correlation";
    string magnetization_folder = comparative_folder + "/magnetization";
    string covariance_folder = comparative_folder + "/covariance";
    string sisj_folder = comparative_folder + "/sisj";
    string sisjsk_folder = comparative_folder + "/sisjsk";
    string triplet_folder = comparative_folder + "/triplet";
    string triplet_same_space_folder = comparative_folder + "/triplet_same_space";
    string triplet_same_number_points_folder = comparative_folder + "/triplet_same_number_points";

    fs::create_directories(correlation_folder);
    fs::create_directories(covariance_folder);
    fs::create_directories(magnetization_folder);
    fs::create_directories(sisj_folder);
    fs::create_directories(sisjsk_folder);
    fs::create_directories(triplet_folder);
    fs::create_directories(triplet_same_space_folder);
    fs::create_directories(triplet_same_number_points_folder);

    // SeparateData subfolders
    string Cij_exp_folder = SeparateData_folder + "/Cij-exp";
    string Cij_ising_folder = SeparateData_folder + "/Cij-ising";
    string h_by_year_folder = SeparateData_folder + "/h_by_year";
    string hi_sp_folder = SeparateData_folder + "/hi";
    string Jij_sp_folder = SeparateData_folder + "/Jij";
    string mi_exp_folder = SeparateData_folder + "/mi-exp";
    string mi_ising_folder = SeparateData_folder + "/mi-ising";
    string Pij_exp_folder = SeparateData_folder + "/Pij-exp";
    string Pij_ising_folder = SeparateData_folder + "/Pij-ising";
    string sisj_exp_folder = SeparateData_folder + "/sisj-exp";
    string sisj_ising_folder = SeparateData_folder + "/sisj-ising";
    string sisjsk_exp_folder = SeparateData_folder + "/sisjsk-exp";
    string sisjsk_ising_folder = SeparateData_folder + "/sisjsk-ising";
    string Tijk_exp_folder = SeparateData_folder + "/Tijk-exp";
    string Tijk_ising_folder = SeparateData_folder + "/Tijk-ising";

    fs::create_directories(Cij_exp_folder);
    fs::create_directories(Cij_ising_folder);
    fs::create_directories(h_by_year_folder);
    fs::create_directories(hi_sp_folder);
    fs::create_directories(Jij_sp_folder);
    fs::create_directories(mi_exp_folder);
    fs::create_directories(mi_ising_folder);
    fs::create_directories(Pij_exp_folder);
    fs::create_directories(Pij_ising_folder);
    fs::create_directories(sisj_exp_folder);
    fs::create_directories(sisj_ising_folder);
    fs::create_directories(sisjsk_exp_folder);
    fs::create_directories(sisjsk_ising_folder);
    fs::create_directories(Tijk_exp_folder);
    fs::create_directories(Tijk_ising_folder);

    // Create Data folders
    string Mag_Corr_folder = "../Data/Mag_Corr";
    fs::create_directories(Mag_Corr_folder);
    
}