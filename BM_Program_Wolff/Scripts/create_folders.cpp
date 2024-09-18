#include "./include/create_folders.h"

int main (){
    // Create folders
	try {
		create_folders();
		cout << "Folders created successfully!" << endl;
    } catch (const fs::filesystem_error& ex) {
        cerr << ex.what() << endl;
    }
}