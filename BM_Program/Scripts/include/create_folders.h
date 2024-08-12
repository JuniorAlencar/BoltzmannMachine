#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;
namespace fs = boost::filesystem;

void create_folders() {
    // Folders exact_solutions-----------------------------------------
    
    // Set folder names using std::string
    string results_folder = "../Results";
    string specificHeat_folder = results_folder + "/SpecificHeat";
    string comparative_folder = results_folder + "/Comparative";
    string CorrJij_folder = results_folder + "/CorrJij";
    string Energy_folder = results_folder + "/Energy";
    string Erro_folder = results_folder + "/Erro";
    string Histogram_folder = results_folder + "/Histogram";
    string Mag_Corr_ising_folder = results_folder + "/Mag_Corr_ising";
    string Magnetization_vs_T_folder = results_folder + "/Magnetization_vs_T";
    string MatrixJij_folder = results_folder + "/MatrixJij";
    string Network_folder = results_folder + "/Network";
    string PJij_folder = results_folder + "/PJij";
    string SeparateData_folder = results_folder + "/SeparateData";

    // Create folders
    fs::create_directories(results_folder);
    fs::create_directories(specificHeat_folder);
    fs::create_directories(comparative_folder);
    fs::create_directories(CorrJij_folder);
    fs::create_directories(Energy_folder);
    fs::create_directories(Erro_folder);
    fs::create_directories(Histogram_folder);
    fs::create_directories(Mag_Corr_ising_folder);
    fs::create_directories(Magnetization_vs_T_folder);
    fs::create_directories(MatrixJij_folder);
    fs::create_directories(Network_folder);
    fs::create_directories(PJij_folder);
    fs::create_directories(SeparateData_folder);

    // Create subfolders
    string hi_folder = Histogram_folder + "/hi";
    string Jij_folder = Histogram_folder + "/Jij";
    string mi_folder = Histogram_folder + "/mi";
    string Pij_folder = Histogram_folder + "/Pij";
    string Tijk_folder = Histogram_folder + "/Tijk";

    fs::create_directories(hi_folder);
    fs::create_directories(Jij_folder);
    fs::create_directories(mi_folder);
    fs::create_directories(Pij_folder);
    fs::create_directories(Tijk_folder);

    // Add more folders as required
    string correlation_folder = comparative_folder + "/correlation";
    string magnetization_folder = comparative_folder + "/magnetization";
    string sisj_folder = comparative_folder + "/sisj";
    string sisjsk_folder = comparative_folder + "/sisjsk";
    string triplet_folder = comparative_folder + "/triplet";
    string triplet_same_space_folder = comparative_folder + "/triplet_same_space";
    string triplet_same_number_points_folder = comparative_folder + "/triplet_same_number_points";

    fs::create_directories(correlation_folder);
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

    // Folders metropolis_solutions-----------------------------------------

    // Set folder names using std::string
    string results_folder_metropolis = "../Results_Metropolis";
    string specificHeat_folder_metropolis = results_folder_metropolis + "/SpecificHeat";
    string comparative_folder_metropolis = results_folder_metropolis + "/Comparative";
    string CorrJij_folder_metropolis = results_folder_metropolis + "/CorrJij";
    string Energy_folder_metropolis = results_folder_metropolis + "/Energy";
    string Erro_folder_metropolis = results_folder_metropolis + "/Erro";
    string Histogram_folder_metropolis = results_folder_metropolis + "/Histogram";
    string Mag_Corr_ising_folder_metropolis = results_folder_metropolis + "/Mag_Corr_ising";
    string Magnetization_vs_T_folder_metropolis = results_folder_metropolis + "/Magnetization_vs_T";
    string MatrixJij_folder_metropolis = results_folder_metropolis + "/MatrixJij";
    string Network_folder_metropolis = results_folder_metropolis + "/Network";
    string PJij_folder_metropolis = results_folder_metropolis + "/PJij";
    string SeparateData_folder_metropolis = results_folder_metropolis + "/SeparateData";

    // Create folders
    fs::create_directories(results_folder_metropolis);
    fs::create_directories(specificHeat_folder_metropolis);
    fs::create_directories(comparative_folder_metropolis);
    fs::create_directories(CorrJij_folder_metropolis);
    fs::create_directories(Energy_folder_metropolis);
    fs::create_directories(Erro_folder_metropolis);
    fs::create_directories(Histogram_folder_metropolis);
    fs::create_directories(Mag_Corr_ising_folder_metropolis);
    fs::create_directories(Magnetization_vs_T_folder_metropolis);
    fs::create_directories(MatrixJij_folder_metropolis);
    fs::create_directories(Network_folder_metropolis);
    fs::create_directories(PJij_folder_metropolis);
    fs::create_directories(SeparateData_folder_metropolis);

    // Create subfolders
    string hi_folder_metropolis = Histogram_folder_metropolis + "/hi";
    string Jij_folder_metropolis = Histogram_folder_metropolis + "/Jij";
    string mi_folder_metropolis = Histogram_folder_metropolis + "/mi";
    string Pij_folder_metropolis = Histogram_folder_metropolis + "/Pij";
    string Tijk_folder_metropolis = Histogram_folder_metropolis + "/Tijk";

    fs::create_directories(hi_folder_metropolis);
    fs::create_directories(Jij_folder_metropolis);
    fs::create_directories(mi_folder_metropolis);
    fs::create_directories(Pij_folder_metropolis);
    fs::create_directories(Tijk_folder_metropolis);

    // Add more folders as required
    string correlation_folder_metropolis = comparative_folder_metropolis + "/correlation";
    string magnetization_folder_metropolis = comparative_folder_metropolis + "/magnetization";
    string sisj_folder_metropolis = comparative_folder_metropolis + "/sisj";
    string sisjsk_folder_metropolis = comparative_folder_metropolis + "/sisjsk";
    string triplet_folder_metropolis = comparative_folder_metropolis + "/triplet";
    string triplet_same_space_folder_metropolis = comparative_folder_metropolis + "/triplet_same_space";
    string triplet_same_number_points_folder_metropolis = comparative_folder_metropolis + "/triplet_same_number_points";

    fs::create_directories(correlation_folder_metropolis);
    fs::create_directories(magnetization_folder_metropolis);
    fs::create_directories(sisj_folder_metropolis);
    fs::create_directories(sisjsk_folder_metropolis);
    fs::create_directories(triplet_folder_metropolis);
    fs::create_directories(triplet_same_space_folder_metropolis);
    fs::create_directories(triplet_same_number_points_folder_metropolis);

    // SeparateData subfolders
    string Cij_exp_folder_metropolis = SeparateData_folder_metropolis + "/Cij-exp";
    string Cij_ising_folder_metropolis = SeparateData_folder_metropolis + "/Cij-ising";
    string h_by_year_folder_metropolis = SeparateData_folder_metropolis + "/h_by_year";
    string hi_sp_folder_metropolis = SeparateData_folder_metropolis + "/hi";
    string Jij_sp_folder_metropolis = SeparateData_folder_metropolis + "/Jij";
    string mi_exp_folder_metropolis = SeparateData_folder_metropolis + "/mi-exp";
    string mi_ising_folder_metropolis = SeparateData_folder_metropolis + "/mi-ising";
    string Pij_exp_folder_metropolis = SeparateData_folder_metropolis + "/Pij-exp";
    string Pij_ising_folder_metropolis = SeparateData_folder_metropolis + "/Pij-ising";
    string sisj_exp_folder_metropolis = SeparateData_folder_metropolis + "/sisj-exp";
    string sisj_ising_folder_metropolis = SeparateData_folder_metropolis + "/sisj-ising";
    string sisjsk_exp_folder_metropolis = SeparateData_folder_metropolis + "/sisjsk-exp";
    string sisjsk_ising_folder_metropolis = SeparateData_folder_metropolis + "/sisjsk-ising";
    string Tijk_exp_folder_metropolis = SeparateData_folder_metropolis + "/Tijk-exp";
    string Tijk_ising_folder_metropolis = SeparateData_folder_metropolis + "/Tijk-ising";

    fs::create_directories(Cij_exp_folder_metropolis);
    fs::create_directories(Cij_ising_folder_metropolis);
    fs::create_directories(h_by_year_folder_metropolis);
    fs::create_directories(hi_sp_folder_metropolis);
    fs::create_directories(Jij_sp_folder_metropolis);
    fs::create_directories(mi_exp_folder_metropolis);
    fs::create_directories(mi_ising_folder_metropolis);
    fs::create_directories(Pij_exp_folder_metropolis);
    fs::create_directories(Pij_ising_folder_metropolis);
    fs::create_directories(sisj_exp_folder_metropolis);
    fs::create_directories(sisj_ising_folder_metropolis);
    fs::create_directories(sisjsk_exp_folder_metropolis);
    fs::create_directories(sisjsk_ising_folder_metropolis);
    fs::create_directories(Tijk_exp_folder_metropolis);
    fs::create_directories(Tijk_ising_folder_metropolis);
}