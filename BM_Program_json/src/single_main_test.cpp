
#include <iostream>
#include "create_folders.h"
#include "exp_means.h"

#include "create_folders.h"

int main(){
	string filename	= "sampleCarmona";
	double min_erro_j	= 1.0e-04;
	double min_erro_h	= 2.0e-04;
	int multiply_teq 	= 300;
	int multiply_relx 	= 2;
	string method = "exact";
	std::cout << "Running crate_Folders test..." << std::endl;	
	c_folders cf;
	int tp = 0;
	string ev;
	ev = cf.create_folders(filename, multiply_teq, multiply_relx, method, tp);
	cout << min_erro_h << endl;
	cout << min_erro_j << endl;
	cout << method << endl;
	std::cout << "Running exp_means test..." << std::endl;
	exp_mean_calculate emc;
	exp_means ep;
	
	ep = emc.exp_calculate(filename);
	return 0;
};



