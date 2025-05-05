#include "./include/create_folders.h"

int main (int argc, char *argv[]){
    string method = argv[1];
    string N_spins = argv[2];
    // Create folders
	try {
		create_folders(method, N_spins);
		cout << "Folders created successfully!" << endl;
    } catch (const fs::filesystem_error& ex) {
        cerr << ex.what() << endl;
    }
}