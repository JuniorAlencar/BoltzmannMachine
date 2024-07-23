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
}