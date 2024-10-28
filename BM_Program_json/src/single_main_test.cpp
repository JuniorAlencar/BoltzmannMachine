
#include <iostream>
#include "create_folders.h"
#include "exp_means.h"
#include "create_folders.h"

int main(int argc, char *argv[]){
	string file_input	= argv[1];
	std::cout << "Running exp_means test..." << std::endl;
	exp_mean_calculate emc;
	exp_means ep;
	
	ep = emc.exp_calculate(file_input);
	return 0;
};



