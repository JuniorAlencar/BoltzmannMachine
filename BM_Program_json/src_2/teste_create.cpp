#include "create_folders.h"

int main(int argc, char *argv[]){
	string filename	= argv[1];
	double min_erro_j	= std::stod(argv[2]);
	double min_erro_h	= std::stod(argv[3]);
	int multiply_teq 	= std::stoi(argv[4]);
	int multiply_relx 	= std::stoi(argv[5]);
	string method = argv[6];
	
	// string filename	= "sampleCarmona";
	// double min_erro_j	= 1.0e-04;
	// double min_erro_h	= 2.0e-04;
	// int multiply_teq 	= 300;
	// int multiply_relx 	= 2;
	// string method = "exact";
	
	c_folders cf;
	    int tp = 0;
	    string ep;
	    ep = cf.create_folders(filename, multiply_teq, multiply_relx, method, tp);
	cout << min_erro_h << endl;
	cout << min_erro_j << endl;
}
