#include "create_folders.h"

int main(){
	string filename	= "sampleCarmona";
	double min_erro_j	= 1.0e-04;
	double min_erro_h	= 2.0e-04;
	int multiply_teq 	= 300;
	int multiply_relx 	= 2;
	string method = "exact";
	c_folders cf;
	    int tp = 0;
	    string ep;
	    ep = cf.create_folders(filename, multiply_teq, multiply_relx, method, tp);
}
