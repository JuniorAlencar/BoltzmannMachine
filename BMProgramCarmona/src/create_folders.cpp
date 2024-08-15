#include "create_folders.h"

#include <boost/filesystem.hpp>
#include <string>

using namespace std;
namespace fs = boost::filesystem;

void create_folders(){
    // Folders to exact_solutions
    string results_folder = "./Results";
    
    fs::create_directories(results_folder);    
    // Folders to metropolis solutions
    string results_metropolis_folder = "./Results_Metropolis";

    
    fs::create_directories(results_metropolis_folder);
}

