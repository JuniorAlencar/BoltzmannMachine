#ifndef NETWORK_H
#define NETWORK_H

#include "nr3.h"

//Classe para criação de uma rede aleatoria. Deve-se entrar com o tamano da rede n, valor medio da distribuição mean, desvio da distribuição sigma e a constante k para a probabilidade de haver ligação
class Rede
{
	public:
		int n, nbonds;
		double mean, sigma; 
		double k;
		VecInt no, s, nb, s_nb;
		VecDoub J;
		VecDoub h;
		
		int type;
		double H;
		
		Rede(int m, double mmean, double ssigma, double kk, int tp, double HH);
		void create_bonds_random (void);
		void create_bonds_tree (void);
		void neighbours(void);
		double gaussian (void);
		void node (int bond);

};

#endif 
