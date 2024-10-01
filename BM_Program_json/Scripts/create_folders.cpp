#include "./include/create_folders.h"

int main (int argc, char *argv[]){
    // Create folders
	try {
        // Arguments to MC
        string text_name	= argv[1];
	    int multiply_teq 	= std::stoi(argv[2]);
	    int multiply_relx 	= std::stoi(argv[3]);
        bool use_exact = (std::string(argv[4]) == "true");
        int type = 0;
		create_folders(text_name, multiply_teq, multiply_relx, use_exact, type);
		cout << "Folders created successfully!" << endl;
    } catch (const fs::filesystem_error& ex) {
        cerr << ex.what() << endl;
    }
}