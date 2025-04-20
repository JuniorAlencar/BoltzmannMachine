#ifndef FORWARDMETHOD_H
#define FORWARDMETHOD_H

#include <vector>
#include <random>
#include <queue>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <climits>

using namespace std;

// Gera um número aleatório entre 0 e 1
inline double random_double() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
}

//CALCULA A ENERGIA TOTAL DO ESTADO ATUAL----------------------------------------------------------
inline double Energy(Rede &r)
{
	double E = 0;
	int k = 0;

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

inline double delta_E(Rede &r, const int s_flip)
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
//-------------------------------------------------------------------------------------------------

inline double Magnetization(Rede &r)
{
	double m = 0;
	int i;
	
	for (i = 0; i < r.n; i++)
	{
		m += r.s[i];
	}
	
	return (double) m;
}

//EXACT SOLUTION-----------------------------------------------------------------------------------
inline void exact_solution(Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const double beta)
{
	int i, n_flip;

	double Z = 0,E0, E_aux, P = 0;
	int state = 0;
	int ind_ss;
	
	VecInt min(r.n, -1.0), max(r.n, -1.0), n_invert(r.n, -1.0);
	int f, alt = 0, b = 1, k, t = 0;

	Z = state = alt = f = 0;

	for (i = 0; i < r.n; i++)
		r.s[i] = 1;		
	
	E0 = Energy(r);
	P = exp(-beta*E0);	
	Z += P;
	
	for (int j = 0; j < r.n; j++)
		av_s[j] += r.s[j]*P;

	ind_ss = 0;
	for (int j = 0; j < r.n; j++)
	{
		for (int l = j; l < r.n; l++)
		{
			av_ss[ind_ss] += r.s[j]*r.s[l]*P;
			ind_ss++;
		}
	}

	for (n_flip = 1; n_flip < r.n; n_flip++)
	{
		for (i = 0; i < r.n; i++)
			r.s[i] = 1;

		for (i = 0; i < n_flip; i++)
		{
			max[i] = r.n - n_flip + i;
			min[i] = i;
		}

		for (i = 0; i < n_flip; i++)
		{
			n_invert[i] = i;
			r.s[i] = -1;
		}

		f = 0;

		while (f == 0)	
		{
			for (i = 0; i < r.n; i++)
				r.s[i] = 1;
			
			E_aux = E0; 

			for (i = 0; i < n_flip; i++)
			{
				E_aux += delta_E (r, n_invert[i]);			
				r.s[n_invert[i]] = -1;
			}
			
			P   = exp(-beta*E_aux);			
			Z  += P;
			
			for (int j = 0; j < r.n; j++)
				av_s[j] += r.s[j]*P;

			ind_ss = 0;
			for (int j = 0; j < r.n; j++)
			{
				for (int l = j; l < r.n; l++)
				{
					av_ss[ind_ss] += r.s[j]*r.s[l]*P;
					ind_ss++;
				}
			}
			
			if (n_invert[0] != max[0])
			{
				b = 1;
				for (i = n_flip-1; i > 0; i--)
				{
					if (n_invert[i] == max[i] && b == 1)
					{				
						if (max[i] == min[i])
						{
							k = 0;
							for (int j = i-1; j > 0; j--)
							{
								if (max[j] != min[j])
								{
									min[i] = min[j] + 2 + k;
									break;
								}
								k++;
							}								
						}
						else 
							min[i]++;
	
						n_invert[i] = min[i];
		
						if (n_invert[i-1] != max[i-1])
						{
							n_invert[i-1]++;
							b = 0;
						}

						alt = 1;
					}
					else
						b = 1;
				}
			
				if (alt == 1)
					alt = 0;
				else
					n_invert[n_flip-1]++;
				
			}
			else
				f = 1;
		}
	}

	for (i = 0; i < r.n; i++)
		r.s[i] = -1.0;

	E_aux = Energy(r);
		
	P = exp(-beta*E_aux);
	Z  += P;
	
	for (int j = 0; j < r.n; j++)
	{
		av_s[j] += r.s[j]*P;
		av_s[j] /= Z;
	}

	ind_ss = 0;
	for (int j = 0; j < r.n; j++)
	{
		for (int l = j; l < r.n; l++)
		{
			av_ss[ind_ss] += r.s[j]*r.s[l]*P;
			av_ss[ind_ss] /= Z;
			ind_ss++;
		}
	}	
}

//EXACT SOLUTION-----------------------------------------------------------------------------------
inline void exact_solution_bm(Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const double beta)
{
	int i, n_flip;

	double Z = 0,E0, E_aux, P = 0;
	int state = 0;
	int ind_ss;
	
	VecInt min(r.n, -1.0), max(r.n, -1.0), n_invert(r.n, -1.0);
	int f, alt = 0, b = 1, k, t = 0;

	Z = state = alt = f = 0;

	for (i = 0; i < r.n; i++)
		r.s[i] = 1;		
	
	E0 = Energy(r);
	P = exp(-beta*E0);	
	Z += P;
	
	for (int j = 0; j < r.n; j++)
		av_s[j] += r.s[j]*P;

	ind_ss = 0;
	for (int j = 0; j < r.n-1; j++)
	{
		for (int l = j+1; l < r.n; l++)
		{
			av_ss[ind_ss] += r.s[j]*r.s[l]*P;
			ind_ss++;
		}
	}

	for (n_flip = 1; n_flip < r.n; n_flip++)
	{
		for (i = 0; i < r.n; i++)
			r.s[i] = 1;

		for (i = 0; i < n_flip; i++)
		{
			max[i] = r.n - n_flip + i;
			min[i] = i;
		}

		for (i = 0; i < n_flip; i++)
		{
			n_invert[i] = i;
			r.s[i] = -1;
		}

		f = 0;

		while (f == 0)	
		{
			for (i = 0; i < r.n; i++)
				r.s[i] = 1;
			
			E_aux = E0; 

			for (i = 0; i < n_flip; i++)
			{
				E_aux += delta_E (r, n_invert[i]);			
				r.s[n_invert[i]] = -1;
			}
			
			P   = exp(-beta*E_aux);			
			Z  += P;
			
			for (int j = 0; j < r.n; j++)
				av_s[j] += r.s[j]*P;

			ind_ss = 0;
			for (int j = 0; j < r.n-1; j++)
			{
				for (int l = j+1; l < r.n; l++)
				{
					av_ss[ind_ss] += r.s[j]*r.s[l]*P;
					ind_ss++;
				}
			}
			
			if (n_invert[0] != max[0])
			{
				b = 1;
				for (i = n_flip-1; i > 0; i--)
				{
					if (n_invert[i] == max[i] && b == 1)
					{				
						if (max[i] == min[i])
						{
							k = 0;
							for (int j = i-1; j > 0; j--)
							{
								if (max[j] != min[j])
								{
									min[i] = min[j] + 2 + k;
									break;
								}
								k++;
							}								
						}
						else 
							min[i]++;
	
						n_invert[i] = min[i];
		
						if (n_invert[i-1] != max[i-1])
						{
							n_invert[i-1]++;
							b = 0;
						}

						alt = 1;
					}
					else
						b = 1;
				}
			
				if (alt == 1)
					alt = 0;
				else
					n_invert[n_flip-1]++;
				
			}
			else
				f = 1;
		}
	}

	for (i = 0; i < r.n; i++)
		r.s[i] = -1.0;

	E_aux = Energy(r);
		
	P = exp(-beta*E_aux);
	Z  += P;
	
	for (int j = 0; j < r.n; j++)
	{
		av_s[j] += r.s[j]*P;
		av_s[j] /= Z;
	}

	ind_ss = 0;
	for (int j = 0; j < r.n-1; j++)
	{
		for (int l = j+1; l < r.n; l++)
		{
			av_ss[ind_ss] += r.s[j]*r.s[l]*P;
			av_ss[ind_ss] /= Z;
			ind_ss++;
		}
	}
}


//EXACT SOLUTION COMPLETA--------------------------------------------------------------------------
inline void exact_solution_comp(Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, Doub &E, Doub &E2, const double beta)
{
	int i, n_flip;

	double Z = 0, E0, E_aux, P = 0;
	int state = 0;
	int ind_ss;
	
	VecInt min(r.n, -1.0), max(r.n, -1.0), n_invert(r.n, -1.0);
	int f, alt = 0, b = 1, k, t = 0;

	Z = state = alt = f = 0;

	for (i = 0; i < r.n; i++)
		r.s[i] = 1;		
	
	E0 = Energy(r);
	P = exp(-beta*E0);	
	
	E += E0*P;
	E2 += E0*E0*P;
	
	Z += P;
	
	for (int j = 0; j < r.n; j++)
		av_s[j] += r.s[j]*P;

	ind_ss = 0;
	for (int j = 0; j < r.n; j++)
	{
		for (int l = j; l < r.n; l++)
		{
			av_ss[ind_ss] += r.s[j]*r.s[l]*P;
			ind_ss++;
		}
	}

	for (n_flip = 1; n_flip < r.n; n_flip++)
	{
		for (i = 0; i < r.n; i++)
			r.s[i] = 1;

		for (i = 0; i < n_flip; i++)
		{
			max[i] = r.n - n_flip + i;
			min[i] = i;
		}

		for (i = 0; i < n_flip; i++)
		{
			n_invert[i] = i;
			r.s[i] = -1;
		}

		f = 0;

		while (f == 0)	
		{
			for (i = 0; i < r.n; i++)
				r.s[i] = 1;
			
			E_aux = E0; 

			for (i = 0; i < n_flip; i++)
			{
				E_aux += delta_E (r, n_invert[i]);			
				r.s[n_invert[i]] = -1;
			}
			
			P   = exp(-beta*E_aux);	
			
			E  += E_aux*P;
			E2 += pow(E_aux, 2)*P;
					
			Z  += P;
			
			for (int j = 0; j < r.n; j++)
				av_s[j] += r.s[j]*P;

			ind_ss = 0;
			for (int j = 0; j < r.n; j++)
			{
				for (int l = j; l < r.n; l++)
				{
					av_ss[ind_ss] += r.s[j]*r.s[l]*P;
					ind_ss++;
				}
			}
			
			if (n_invert[0] != max[0])
			{
				b = 1;
				for (i = n_flip-1; i > 0; i--)
				{
					if (n_invert[i] == max[i] && b == 1)
					{				
						if (max[i] == min[i])
						{
							k = 0;
							for (int j = i-1; j > 0; j--)
							{
								if (max[j] != min[j])
								{
									min[i] = min[j] + 2 + k;
									break;
								}
								k++;
							}								
						}
						else 
							min[i]++;
	
						n_invert[i] = min[i];
		
						if (n_invert[i-1] != max[i-1])
						{
							n_invert[i-1]++;
							b = 0;
						}

						alt = 1;
					}
					else
						b = 1;
				}
			
				if (alt == 1)
					alt = 0;
				else
					n_invert[n_flip-1]++;
				
			}
			else
				f = 1;
		}
	}

	for (i = 0; i < r.n; i++)
		r.s[i] = -1.0;

	E_aux = Energy(r);
		
	P = exp(-beta*E_aux);
	
	E  += E_aux*P;
	E2 += pow(E_aux, 2)*P;
	
	Z  += P;
	
	for (int j = 0; j < r.n; j++)
	{
		av_s[j] += r.s[j]*P;
		av_s[j] /= Z;
	}

	ind_ss = 0;
	for (int j = 0; j < r.n; j++)
	{
		for (int l = j; l < r.n; l++)
		{
			av_ss[ind_ss] += r.s[j]*r.s[l]*P;
			av_ss[ind_ss] /= Z;
			ind_ss++;
		}
	}
	
	E  /= Z;
	E2 /= Z;
	
}

//METROPOLIS---------------------------------------------------------------------------------------
inline void metropolis(Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const int t_eq, const int t_step,  
				const int relx, const int rept, const double beta)
{
	int s_flip;
	double p;
	double dE;
	int ind_ss = 0;
	
	
	for (int j = 0; j < rept; j++)
	{
		
		for (int i = 0; i < t_eq + t_step; i++)
		{
			if (i >= t_eq && (int) (i-t_eq)%(int)relx == 0)
			{
				for (int jj = 0; jj < r.n; jj++)
					av_s[jj] += r.s[jj];

				ind_ss = 0;
				for (int jj = 0; jj < r.n; jj++)
				{
					for (int l = jj; l < r.n; l++)
					{
						av_ss[ind_ss] += r.s[jj]*r.s[l];
						ind_ss++;
					}
				}
			}
	
			s_flip = rand() % r.n; 
			dE = delta_E (r, s_flip);
			
			if (dE <= 0)
			{
				r.s[s_flip] = -r.s[s_flip];
			}
			else
			{
				 p = ((double) rand()/RAND_MAX);
				 
				 if (p < exp(-beta*dE))
				 {
				 	r.s[s_flip] = -r.s[s_flip];
				 }
			}
		}
	
	}
	
	for (int jj = 0; jj < r.n; jj++)
		av_s[jj] /= (rept*t_step/relx);

	ind_ss = 0;
	for (int jj = 0; jj < r.n; jj++)
	{
		for (int l = jj; l < r.n; l++)
		{
			av_ss[ind_ss] /= (rept*t_step/relx);
			ind_ss++;
		}
	}

}

//METROPOLIS PARA BM-------------------------------------------------------------------------------
inline void metropolis_bm(Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const int t_eq, const int t_step,  
					const int relx, const int rept, const double beta)
{
	int s_flip;
	double p;
	double dE;
	int ind_ss = 0;
	
	
	for (int j = 0; j < rept; j++)
	{
		
		for (int i = 0; i < t_eq + t_step; i++)
		{
			if (i >= t_eq && (int) (i-t_eq)%(int)relx == 0)
			{
				for (int jj = 0; jj < r.n; jj++)
					av_s[jj] += r.s[jj];

				ind_ss = 0;
				for (int jj = 0; jj < r.n-1; jj++)
				{
					for (int l = jj+1; l < r.n; l++)
					{
						av_ss[ind_ss] += r.s[jj]*r.s[l];
						ind_ss++;
					}
				}
			}
	
			s_flip = rand() % r.n; 
			dE = delta_E (r, s_flip);
			
			if (dE <= 0)
			{
				r.s[s_flip] = -r.s[s_flip];
			}
			else
			{
				 p = ((double) rand()/RAND_MAX);
				 
				 if (p < exp(-beta*dE))
				 {
				 	r.s[s_flip] = -r.s[s_flip];
				 }
			}
		}
	
	}
	
	for (int jj = 0; jj < r.n; jj++)
		av_s[jj] /= (rept*t_step/relx);

	ind_ss = 0;
	for (int jj = 0; jj < r.n-1; jj++)
	{
		for (int l = jj+1; l < r.n; l++)
		{
			av_ss[ind_ss] /= (rept*t_step/relx);
			ind_ss++;
		}
	}
}

//METROPOLIS COMPLETO------------------------------------------------------------------------------
inline void metropolis_comp(Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const int t_eq, const int t_step,  const int relx, const int rept, const double beta, Doub &E, Doub &E2)
{
	int s_flip;
	double p;
	double dE;
	double E0;
	int ind_ss = 0;
	
	E0 = Energy(r);
	
	for (int j = 0; j < rept; j++)
	{
		
		for (int i = 0; i < t_eq + t_step; i++)
		{
			if (i >= t_eq && (int) (i-t_eq)%(int)relx == 0)
			{
				for (int jj = 0; jj < r.n; jj++)
					av_s[jj] += r.s[jj];

				ind_ss = 0;
				for (int jj = 0; jj < r.n; jj++)
				{
					for (int l = jj; l < r.n; l++)
					{
						av_ss[ind_ss] += r.s[jj]*r.s[l];
						ind_ss++;
					}
				}
				
				E  += E0;
				E2 += pow(E0, 2);
			}
	
			s_flip = rand() % r.n; 
			dE = delta_E (r, s_flip);
			
			if (dE <= 0)
			{
				r.s[s_flip] = -r.s[s_flip];
				E0 += dE;
			}
			else
			{
				 p = ((double) rand()/RAND_MAX);
				 
				 if (p < exp(-beta*dE))
				 {
				 	r.s[s_flip] = -r.s[s_flip];
				 	E0 += dE;
				 }
			}
		}
	
	}
	
	for (int jj = 0; jj < r.n; jj++)
		av_s[jj] /= (rept*t_step/relx);

	ind_ss = 0;
	for (int jj = 0; jj < r.n; jj++)
	{
		for (int l = jj; l < r.n; l++)
		{
			av_ss[ind_ss] /= (rept*t_step/relx);
			ind_ss++;
		}
	}
	
	E /= (double) (rept*t_step/relx);
	E2 /= (double) (rept*t_step/relx);
	
}


//-----------------------------------------------------------------------------
//Solução exata para o calor especifico

inline void exact_solution_cap(Rede &r, Doub &E, Doub &E2, const double beta)
{
	int i, n_flip;

	double Z = 0, E0, E_aux, P = 0;
	int state = 0;
	int ind_ss;
	
	VecInt min(r.n, -1.0), max(r.n, -1.0), n_invert(r.n, -1.0);
	int f, alt = 0, b = 1, k, t = 0;

	Z = state = alt = f = 0;

	for (i = 0; i < r.n; i++)
		r.s[i] = 1;		
	
	E0 = Energy(r);
	P = exp(-beta*E0);	
	
	E += E0*P;
	E2 += E0*E0*P;
	
	Z += P;

	for (n_flip = 1; n_flip < r.n; n_flip++)
	{
		for (i = 0; i < r.n; i++)
			r.s[i] = 1;

		for (i = 0; i < n_flip; i++)
		{
			max[i] = r.n - n_flip + i;
			min[i] = i;
		}

		for (i = 0; i < n_flip; i++)
		{
			n_invert[i] = i;
			r.s[i] = -1;
		}

		f = 0;

		while (f == 0)	
		{
			for (i = 0; i < r.n; i++)
				r.s[i] = 1;
			
			E_aux = E0; 

			for (i = 0; i < n_flip; i++)
			{
				E_aux += delta_E (r, n_invert[i]);			
				r.s[n_invert[i]] = -1;
			}
			
			P   = exp(-beta*E_aux);	
			
			E  += E_aux*P;
			E2 += pow(E_aux, 2)*P;
					
			Z  += P;
			
			if (n_invert[0] != max[0])
			{
				b = 1;
				for (i = n_flip-1; i > 0; i--)
				{
					if (n_invert[i] == max[i] && b == 1)
					{				
						if (max[i] == min[i])
						{
							k = 0;
							for (int j = i-1; j > 0; j--)
							{
								if (max[j] != min[j])
								{
									min[i] = min[j] + 2 + k;
									break;
								}
								k++;
							}								
						}
						else 
							min[i]++;
	
						n_invert[i] = min[i];
		
						if (n_invert[i-1] != max[i-1])
						{
							n_invert[i-1]++;
							b = 0;
						}

						alt = 1;
					}
					else
						b = 1;
				}
			
				if (alt == 1)
					alt = 0;
				else
					n_invert[n_flip-1]++;
				
			}
			else
				f = 1;
		}
	}

	for (i = 0; i < r.n; i++)
		r.s[i] = -1.0;

	E_aux = Energy(r);
		
	P = exp(-beta*E_aux);
	
	E  += E_aux*P;
	E2 += pow(E_aux, 2)*P;
	
	Z  += P;
	
	E  /= Z;
	E2 /= Z;
	
}

//METROPOLIS para o calor especifico-------------------------------------------
inline void metropolis_cap(Rede &r, const int t_eq, const int t_step,  const int relx, const int rept, const double T, 
						Doub &E, Doub &E2)
{
	int s_flip;
	double p;
	double dE;
	double E0;
	
	E0 = Energy(r);
	
	for (int j = 0; j < rept; j++)
	{
		
		for (int i = 0; i < t_eq + t_step; i++)
		{
			if (i >= t_eq && (int) (i-t_eq)%(int)relx == 0)
			{				
				E  += E0;
				E2 += pow(E0, 2);
			}
	
			s_flip = rand() % r.n; 
			dE = delta_E (r, s_flip);
			
			if (dE <= 0)
			{
				r.s[s_flip] = -r.s[s_flip];
				E0 += dE;
			}
			else
			{
				 p = ((double) rand()/RAND_MAX);
				 
				 if (p < exp(-dE/T))
				 {
				 	r.s[s_flip] = -r.s[s_flip];
				 	E0 += dE;
				 }
			}
		}
	
	}
	
	E /= (rept*t_step/relx);
	E2 /= (rept*t_step/relx);
	
}


//METROPOLIS para o calor especifico e magnetização-------------------------------------------
inline void metropolis_cap_mag(Rede &r, const int t_eq, const int t_step,  const int relx, const int rept, const double T, 
						Doub &E, Doub &E2, double &M, double &M2)
{
	int s_flip;
	double p;
	double dE;
	double E0;
	double M0;
	
	E0 = Energy(r);
	
	for (int j = 0; j < rept; j++)
	{
		
		for (int i = 0; i < t_eq + t_step; i++)
		{
			if (i >= t_eq && (int) (i-t_eq)%(int)relx == 0)
			{				
				E  += E0;
				E2 += pow(E0, 2);

				M0 = Magnetization(r);
                M += M0;
                M2 += M0*M0;
			}
	
			s_flip = rand() % r.n; 
			dE = delta_E (r, s_flip);
			
			if (dE <= 0)
			{
				r.s[s_flip] = -r.s[s_flip];
				E0 += dE;
			}
			else
			{
				 p = ((double) rand()/RAND_MAX);
				 
				 if (p < exp(-dE/T))
				 {
				 	r.s[s_flip] = -r.s[s_flip];
				 	E0 += dE;
				 }
			}
		}
	
	}
	
	E /= (rept*t_step/relx);
	E2 /= (rept*t_step/relx);

	M /= (rept*t_step/relx);
    M2 /= (rept*t_step/relx);
	
}



//-------------------------------------------------------------------
// Terceiro momentos

//METROPOLIS---------------------------------------------------------------------------------------
inline void metropolis_triplet(Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, vector<double> &av_sss, const int t_eq, const int t_step,  
				const int relx, const int rept, const double beta)
{
	int s_flip;
	double p;
	double dE;
	int ind_ss = 0;
	
	
	for (int j = 0; j < rept; j++)
	{
		
		for (int i = 0; i < t_eq + t_step; i++)
		{
			if (i >= t_eq && (int) (i-t_eq)%(int)relx == 0)
			{
				for (int jj = 0; jj < r.n; jj++)
					av_s[jj] += r.s[jj];

				ind_ss = 0;
				for (int jj = 0; jj < r.n-1; jj++)
				{
					for (int l = jj+1; l < r.n; l++)
					{
						av_ss[ind_ss] += r.s[jj]*r.s[l];
						ind_ss++;
					}
				}

				ind_ss = 0;
				for (int l = 0; l < r.n-2; l++)
				{
					for (int m = l+1; m < r.n-1; m++)
					{
						for (int p = m+1; p < r.n; p++)
						{
							av_sss[ind_ss] += r.s[l]*r.s[m]*r.s[p];
							ind_ss++;
						}
					}
				}
			}
	
			s_flip = rand() % r.n; 
			dE = delta_E (r, s_flip);
			
			if (dE <= 0)
			{
				r.s[s_flip] = -r.s[s_flip];
			}
			else
			{
				 p = ((double) rand()/RAND_MAX);
				 
				 if (p < exp(-beta*dE))
				 {
				 	r.s[s_flip] = -r.s[s_flip];
				 }
			}
		}
	
	}
	
	for (int jj = 0; jj < r.n; jj++)
		av_s[jj] /= (rept*t_step/relx);

	ind_ss = 0;
	for (int jj = 0; jj < r.n-1; jj++)
	{
		for (int l = jj+1; l < r.n; l++)
		{
			av_ss[ind_ss] /= (rept*t_step/relx);
			ind_ss++;
		}
	}

	ind_ss = 0;
	for (int l = 0; l < r.n-2; l++)
	{
		for (int m = l+1; m < r.n-1; m++)
		{
			for (int p = m+1; p < r.n; p++)
			{
				av_sss[ind_ss] /= (rept*t_step/relx);
				ind_ss++;
			}
		}
	}

}


inline void exact_solution_triplet(Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, vector<double> &av_sss, const double beta)
{
	int i, n_flip;

	double Z = 0,E0, E_aux, P = 0;
	int state = 0;
	int ind_ss;
	
	VecInt min(r.n, -1.0), max(r.n, -1.0), n_invert(r.n, -1.0);
	int f, alt = 0, b = 1, k, t = 0;

	Z = state = alt = f = 0;

	for (i = 0; i < r.n; i++)
		r.s[i] = 1;		
	
	E0 = Energy(r);
	P = exp(-beta*E0);	
	Z += P;
	
	for (int j = 0; j < r.n; j++)
		av_s[j] += r.s[j]*P;

	ind_ss = 0;
	for (int j = 0; j < r.n-1; j++)
	{
		for (int l = j+1; l < r.n; l++)
		{
			av_ss[ind_ss] += r.s[j]*r.s[l]*P;
			ind_ss++;

		}
	}

	ind_ss = 0;
	for (int j = 0; j < r.n-2; j++)
	{
		for (int l = j+1; l < r.n-1; l++)
		{
			for (int k = l+1; k < r.n; k++)
			{  
				av_sss[ind_ss] += r.s[j]*r.s[l]*r.s[k]*P;
				ind_ss++;
			}
		}
	}

	for (n_flip = 1; n_flip < r.n; n_flip++)
	{
		for (i = 0; i < r.n; i++)
			r.s[i] = 1;

		for (i = 0; i < n_flip; i++)
		{
			max[i] = r.n - n_flip + i;
			min[i] = i;
		}

		for (i = 0; i < n_flip; i++)
		{
			n_invert[i] = i;
			r.s[i] = -1;
		}

		f = 0;

		while (f == 0)	
		{
			for (i = 0; i < r.n; i++)
				r.s[i] = 1;
			
			E_aux = E0; 

			for (i = 0; i < n_flip; i++)
			{
				E_aux += delta_E (r, n_invert[i]);			
				r.s[n_invert[i]] = -1;
			}
			
			P   = exp(-beta*E_aux);			
			Z  += P;
			
			for (int j = 0; j < r.n; j++)
				av_s[j] += r.s[j]*P;

			ind_ss = 0;
			for (int j = 0; j < r.n-1; j++)
			{
				for (int l = j+1; l < r.n; l++)
				{
					av_ss[ind_ss] += r.s[j]*r.s[l]*P;
					ind_ss++;
				}
			}
			
			ind_ss = 0;
			for (int j = 0; j < r.n-2; j++)
			{
				for (int l = j+1; l < r.n-1; l++)
				{
					for (int k = l+1; k < r.n; k++)
					{  
						av_sss[ind_ss] += r.s[j]*r.s[l]*r.s[k]*P;
						ind_ss++;
					}
				}
			}
			
			if (n_invert[0] != max[0])
			{
				b = 1;
				for (i = n_flip-1; i > 0; i--)
				{
					if (n_invert[i] == max[i] && b == 1)
					{				
						if (max[i] == min[i])
						{
							k = 0;
							for (int j = i-1; j > 0; j--)
							{
								if (max[j] != min[j])
								{
									min[i] = min[j] + 2 + k;
									break;
								}
								k++;
							}								
						}
						else 
							min[i]++;
	
						n_invert[i] = min[i];
		
						if (n_invert[i-1] != max[i-1])
						{
							n_invert[i-1]++;
							b = 0;
						}

						alt = 1;
					}
					else
						b = 1;
				}
			
				if (alt == 1)
					alt = 0;
				else
					n_invert[n_flip-1]++;
				
			}
			else
				f = 1;
		}
	}

	for (i = 0; i < r.n; i++)
		r.s[i] = -1.0;

	E_aux = Energy(r);
		
	P = exp(-beta*E_aux);
	Z  += P;
	
	for (int j = 0; j < r.n; j++)
	{
		av_s[j] += r.s[j]*P;
		av_s[j] /= Z;
	}

	ind_ss = 0;
	for (int j = 0; j < r.n-1; j++)
	{
		for (int l = j+1; l < r.n; l++)
		{
			av_ss[ind_ss] += r.s[j]*r.s[l]*P;
			av_ss[ind_ss] /= Z;
			ind_ss++;
		}
	}
	
	ind_ss = 0;
	for (int j = 0; j < r.n-2; j++)
	{
		for (int l = j+1; l < r.n-1; l++)
		{
			for (int k = l+1; k < r.n; k++)
			{  
				av_sss[ind_ss] += r.s[j]*r.s[l]*r.s[k]*P;
				av_sss[ind_ss] /= Z;
				ind_ss++;
			}
		}
	}
	
}

inline void parallel_tempering(
    int n_replicas, double T_min, double T_max,
    VecDoub_IO &av_s, VecDoub_IO &av_ss,
    int t_eq, int t_step, int relx, int rept,
    int n_spins, double mean, double sigma,
    const int &type, double H,
    std::vector<double> &energy_per_replica,
    std::vector<double> &temperatures,
    double &swap_acceptance_ratio
) {
    // Construir vetor de temperaturas (log-scale)
    temperatures.resize(n_replicas);
    for (int i = 0; i < n_replicas; ++i)
        temperatures[i] = T_min * pow(T_max / T_min, (double)i / (n_replicas - 1));


    // Inicializar réplicas
    std::vector<Rede> replicas;
    for (int i = 0; i < n_replicas; ++i) {
        double beta = 1.0 / temperatures[i];
        Rede r(n_spins, mean, sigma, beta, type, H);
        r.create_bonds_random();
        replicas.push_back(r);
    }

    // Encontrar a réplica mais próxima de β = 1.0
    int target_index = -1;
    double min_diff = 1e9;
    for (int i = 0; i < replicas.size(); ++i) {
        double diff = fabs(replicas[i].k - 1.0);
        if (diff < min_diff) {
            min_diff = diff;
            target_index = i;
        }
    }

    if (target_index == -1) {
        std::cerr << "Nenhuma réplica com beta ≈ 1.0 encontrada." << std::endl;
        return;
    }

    // Zerando acumuladores
    energy_per_replica.assign(n_replicas, 0.0);

    int n_meas = 0;
    int swap_attempts = 0, swap_accepted = 0;

    for (int rep = 0; rep < rept; ++rep) {
        // Equilibração
        for (int step = 0; step < t_eq; ++step) {
            for (auto &replica : replicas) {
                int s_flip = rand() % replica.n;
                double dE = delta_E(replica, s_flip);
                if (dE <= 0 || ((double)rand() / RAND_MAX) < exp(-replica.k * dE))
                    replica.s[s_flip] *= -1;
            }
        }

        // Amostragem
        for (int step = 0; step < t_step; ++step) {
            for (auto &replica : replicas) {
                int s_flip = rand() % replica.n;
                double dE = delta_E(replica, s_flip);
                if (dE <= 0 || ((double)rand() / RAND_MAX) < exp(-replica.k * dE))
                    replica.s[s_flip] *= -1;
            }

            // Acumular médias da réplica-alvo (β ≈ 1.0)
            if (step % relx == 0) {
                Rede &replica = replicas[target_index];
                for (int i = 0; i < replica.n; ++i)
                    av_s[i] += replica.s[i];

                int ind_ss = 0;
                for (int i = 0; i < replica.n - 1; ++i)
                    for (int j = i + 1; j < replica.n; ++j)
                        av_ss[ind_ss++] += replica.s[i] * replica.s[j];

                ++n_meas;
            }

            // Trocas entre réplicas adjacentes
            for (int i = 0; i < n_replicas - 1; ++i) {
                double beta_i = replicas[i].k;
                double beta_j = replicas[i + 1].k;
                double E_i = Energy(replicas[i]);
                double E_j = Energy(replicas[i + 1]);

                double delta = (beta_j - beta_i) * (E_j - E_i);
                ++swap_attempts;
                if ((double)rand() / RAND_MAX < exp(delta)) {
                    std::swap(replicas[i].s, replicas[i + 1].s);
                    ++swap_accepted;
                }

                energy_per_replica[i] += E_i;
                if (i == n_replicas - 2)
                    energy_per_replica[i + 1] += E_j;
            }
        }
    }

    // Verifica se houve medidas antes de normalizar
    if (n_meas > 0) {
        for (int i = 0; i < av_s.size(); ++i)
            av_s[i] /= n_meas;

        int ind_ss = 0;
        for (int i = 0; i < n_spins - 1; ++i)
            for (int j = i + 1; j < n_spins; ++j)
                av_ss[ind_ss++] /= n_meas;
    } else {
        std::cerr << "Aviso: nenhuma medida foi acumulada (n_meas == 0). Verifique a escolha de beta_target." << std::endl;
    }

    for (int i = 0; i < n_replicas; ++i)
        energy_per_replica[i] /= (rept * t_step);

    swap_acceptance_ratio = (swap_attempts > 0) ? (double)swap_accepted / swap_attempts : 0.0;
}

void parallel_tempering_cuda_multi(
    int n_replicas, double T_min, double T_max,
    VecDoub_IO &bm_av_s, VecDoub_IO &bm_av_ss,
    int t_eq, int t_step, int relx, int rept,
    int n_spins, double mean, double sigma,
    const int &type, double H,
    std::vector<double> &energy_per_replica,
    std::vector<double> &temperatures,
    double &swap_acceptance_ratio,
	std::mt19937 &gen
);



#endif
