#ifndef NETWORK_H
#define NETWORK_H

#include "nr3.h"

//Classe para criação de uma rede aleatória. 
// Deve-se entrar com o tamanho da rede n, 
// o valor médio da distribuição mean, 
// o desvio da distribuição sigma 
// e a constante k para a probabilidade de haver ligação
class Rede
{
	public:
		// Number spins and 
		int n, nbonds;
		// 
		double mean, sigma; 
		// Probability of connection
		double k;
		VecInt no, s, nb, s_nb;
		// Vector that will carry the values ​​of Jij
		VecDoub J;
		// Vector that will carry the values ​​of hi
		VecDoub h;
		// Type 0 -> random network, Type 1 -> regular network
		int type;
		double H;
		// Constructor to network
		Rede(const int &m, 
		const double &mean_, 
		const double &sigma_, 
		const double &k_, 
		const int &type_, 
		const double &H_);
		void create_bonds_random (void);
		void create_bonds_tree (void);
		void neighbors(void);
		double gaussian (void);
		void node (int bond);

};


#endif 
