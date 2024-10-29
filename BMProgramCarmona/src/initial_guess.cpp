#include "initial_guess.h"
#include <iostream>
#include <fstream>
using namespace std;

void read_initial_guess(const string filename, Rede& bm){
    int nspins;
    string line;

	ifstream fin(filename.c_str());

	if (!fin) {
    std::cerr << "Unable to open file " << filename << endl;
    exit(1); // or handle the error in a different way
	}


	fin >> nspins;

	if (nspins != bm.n){
		std::cerr << "file not compatible " << filename << endl;
    	exit(1); // or handle the error in a different way

	}
	cout << "lendo h0 e J0 de " << filename << endl;

    int nbonds=nspins*(nspins-1)/2;
	for (int i = 0; i < nbonds; i++)
	{
		fin >> bm.J[i];

		if (i < nspins)
			fin >> bm.h[i];
	}
	fin.close();
}