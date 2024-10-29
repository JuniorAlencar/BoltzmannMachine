#ifndef CREATE_FOLDERS_H
#define CREATE_FOLDERS_H

#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <fmt/core.h>
#include <filesystem>
#include <algorithm>


using namespace std;
namespace fs = std::filesystem;

// variable type is an auxiliary variable, where
// type=0 returns nothing,
// type=1 returns properties folder
// type=2 returns errors folder
// type=3 returns network folder
// type=4 returns specific heat folder

class c_folders{
    public:
        std::string create_folders(const string &text_name, const int &multiply_teq, const int &multiply_relx, const string &method, const int &type);
};

#endif // CREATE_FOLDERS_H
