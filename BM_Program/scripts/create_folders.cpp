#include "./include/create_folders.h"

int main (int argc, char *argv[]){
    string method = argv[1];
    // Create folders
	try {
		create_folders(method);
		cout << "Folders created successfully!" << endl;
    } catch (const fs::filesystem_error& ex) {
        cerr << ex.what() << endl;
    }
}