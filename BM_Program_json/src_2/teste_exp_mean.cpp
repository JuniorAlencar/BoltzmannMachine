#include "exp_means.h"

using namespace std;

int main(int argc, char *argv[]){
	string file_input	= argv[1];
	cout << "Running exp_means test..." << endl;
	exp_mean_calculate emc;
	exp_means ep;
	
	ep = emc.exp_calculate(file_input);
	return 0;
};
