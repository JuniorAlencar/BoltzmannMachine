#include "energy.h"

// Network arguments in network.h

//CALCULA A ENERGIA TOTAL DO ESTADO ATUAL----------------------------------------------------------
double Energy (Rede &r)
{
	// Initial Energy
	double E = 0;
	int k = 0;
	
	// Sum of the contribution of all spins
	for (int i = 0; i < r.nbonds; i++)
	{
		if (i < r.n)
			
			E += r.h[i]*r.s[i];
		
		r.Rede::node(i);
	
		E += r.J[i]*r.s[r.no[0]]*r.s[r.no[1]];		
		
	}

	return -E;
	
}
//-------------------------------------------------------------------------------------------------

// Energy difference after flip spin
double delta_E (Rede &r, const int s_flip)
{
	double dE, Jsum = 0;
	int i;
	
	for (i = s_flip*(r.n-1); i < (s_flip+1)*(r.n-1); i++)
	{
		Jsum += r.J[r.nb[i]]*r.s[r.s_nb[i]];
	}
		
	dE = 2*r.s[s_flip]*(Jsum + r.h[s_flip]);	
	
	return dE;	
}


double Magnetization (Rede &r)
{
	double m = 0;
	int i;
	
	for (i = 0; i < r.n; i++)
	{
		m += r.s[i];
	}
	
	return (double) m;
}