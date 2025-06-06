#ifndef FORWARDMETHOD_H
#define FORWARDMETHOD_H

#include <vector>
#include <random>
#include <queue>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <climits>

using namespace std;

// Gera um número aleatório entre 0 e 1
double random_double() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
}

//CALCULA A ENERGIA TOTAL DO ESTADO ATUAL----------------------------------------------------------
double Energy (Rede &r)
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
//-------------------------------------------------------------------------------------------------

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

//EXACT SOLUTION-----------------------------------------------------------------------------------
void exact_solution (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const double beta)
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
void exact_solution_bm (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const double beta)
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
void exact_solution_comp (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, Doub &E, Doub &E2, const double beta)
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
void metropolis (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const int t_eq, const int t_step,  
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

int draw_site(Rede &r, std::mt19937 &gen) {
    static std::uniform_int_distribution<int> dist(0, r.n - 1);
    return dist(gen);
}

double draw_probability(std::mt19937 &gen) {
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(gen);
}


//METROPOLIS PARA BM-------------------------------------------------------------------------------
void metropolis_bm(Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const int t_eq, const int t_step,  
	const int relx, const int rept, const double beta, std::mt19937 &gen)
{
	int s_flip;
	double p;
	double dE;
	int ind_ss = 0;

	for (int j = 0; j < rept; j++){
			for (int i = 0; i < t_eq + t_step; i++)
		{
			if (i >= t_eq && (i - t_eq) % relx == 0)
		{
			for (int jj = 0; jj < r.n; jj++)
			av_s[jj] += r.s[jj];

		ind_ss = 0;
		for (int jj = 0; jj < r.n - 1; jj++)
		{
			for (int l = jj + 1; l < r.n; l++)
			{
				av_ss[ind_ss] += r.s[jj] * r.s[l];
				ind_ss++;
		}
		}
	}

		s_flip = draw_site(r, gen);
		dE = delta_E(r, s_flip);

		if (dE <= 0){
			r.s[s_flip] = -r.s[s_flip];
		}
		else{
			p = draw_probability(gen);
			if (p < exp(-beta * dE)){
			r.s[s_flip] = -r.s[s_flip];
				}
			}
		}
	}

	for (int jj = 0; jj < r.n; jj++)
		av_s[jj] /= (rept * t_step / relx);

	ind_ss = 0;
	for (int jj = 0; jj < r.n - 1; jj++){
		for (int l = jj + 1; l < r.n; l++){
		av_ss[ind_ss] /= (rept * t_step / relx);
		ind_ss++;
		}
	}
}


//METROPOLIS COMPLETO------------------------------------------------------------------------------
void metropolis_comp (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const int t_eq, const int t_step,  const int relx, const int rept, const double beta, Doub &E, Doub &E2)
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

void exact_solution_cap (Rede &r, Doub &E, Doub &E2, const double beta)
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
void metropolis_cap (Rede &r, const int t_eq, const int t_step,  const int relx, const int rept, const double T, 
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
void metropolis_cap_mag (Rede &r, const int t_eq, const int t_step,  const int relx, const int rept, const double T, 
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
void metropolis_triplet (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, vector<double> &av_sss, const int t_eq, const int t_step,  
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


void exact_solution_triplet (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, vector<double> &av_sss, const double beta)
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

// Adaptative temperatures by Kone & Kofth
// void generate_adaptive_temperatures(double &T_min, double &T_max ,std::vector<double> &temperatures, const std::vector<double> &energies, int n_replicas) {
//     double mean_E = 0.0, var_E = 0.0;
//     int N = energies.size();
//     for (double E : energies) mean_E += E;
//     mean_E /= N;
//     for (double E : energies) var_E += (E - mean_E) * (E - mean_E);
//     var_E /= N;
//     double sigma_E = sqrt(var_E);
	
//     // Ajustar T_min e T_max baseado no sigma_E
//     if (sigma_E > 1e-8) {
//         T_min = 1.0 / (2.0 * sigma_E);  // Mais físico: inversamente proporcional à flutuação
//     } else {
//         T_min = 2.0;  // fallback padrão
//     }
//     T_max = T_min * 7.0;  // ou 10.0x, depende da sua preferência

//     // Agora sim, gerar distribuição logarítmica entre novo T_min e T_max
//     temperatures.resize(n_replicas);
//     for (int i = 0; i < n_replicas; ++i) {
//         double alpha = (double)i / (n_replicas - 1);
//         temperatures[i] = T_min * pow(T_max / T_min, alpha);
//     }
// }

void generate_temperatures_optimized(
    double T_min, double T_max, 
    std::vector<double> &temperatures, 
    int n_replicas,
    double alpha_mix = 0.8 // 0 = linear, 1 = log, intermediário é misto
) {
    temperatures.resize(n_replicas);

    for (int i = 0; i < n_replicas; ++i) {
        double frac = (double)i / (n_replicas - 1);

        // Combinação linear + logarítmica
        double T_log = T_min * pow(T_max / T_min, frac);
        double T_lin = T_min + frac * (T_max - T_min);

        temperatures[i] = alpha_mix * T_log + (1.0 - alpha_mix) * T_lin;
    }
}

void parallel_tempering(
    std::vector<Rede> &replicas,
    VecDoub_IO &av_s, VecDoub_IO &av_ss,
    const int t_eq, const int t_step, const int relx, const int rept,
    int n_spins, double mean, double sigma,
    const int &type, double H,
    std::vector<double> &energy_per_replica,
    double &swap_acceptance_ratio,
    std::mt19937 &gen_global
) {
    int n_replicas = replicas.size();
    int n_meas = 0;
    int swap_attempts = 0, swap_accepted = 0;

    // --- Encontrar réplica alvo (beta mais próximo de 1) ---
    int target_index = 0;
    for (int i = 1; i < n_replicas; ++i)
        if (fabs(replicas[i].k - 1.0) < fabs(replicas[target_index].k - 1.0))
            target_index = i;

    // --- Reset Acumuladores ---
    energy_per_replica.assign(n_replicas, 0.0);

    int num_threads = omp_get_num_procs();
    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
        std::mt19937 gen_local(gen_global() + omp_get_thread_num());

        #pragma omp single
        {
            for (int i = 0; i < n_replicas; ++i) {
                #pragma omp task firstprivate(i)
                {
                    for (int rep = 0; rep < rept; ++rep) {
                        // --- Equilibração ---
                        for (int step = 0; step < t_eq; ++step) {
                            int s_flip = draw_site(replicas[i], gen_local);
                            double dE = delta_E(replicas[i], s_flip);
                            if (dE <= 0 || draw_probability(gen_local) < exp(-replicas[i].k * dE))
                                replicas[i].s[s_flip] *= -1;
                        }

                        // --- Amostragem ---
                        for (int step = 0; step < t_step; ++step) {
                            int s_flip = draw_site(replicas[i], gen_local);
                            double dE = delta_E(replicas[i], s_flip);
                            if (dE <= 0 || draw_probability(gen_local) < exp(-replicas[i].k * dE))
                                replicas[i].s[s_flip] *= -1;

                            if (i == target_index && step % relx == 0) {
                                #pragma omp critical
                                {
                                    for (int j = 0; j < replicas[i].n; ++j)
                                        av_s[j] += replicas[i].s[j];

                                    int ind_ss = 0;
                                    for (int j = 0; j < replicas[i].n - 1; ++j)
                                        for (int k = j + 1; k < replicas[i].n; ++k)
                                            av_ss[ind_ss++] += replicas[i].s[j] * replicas[i].s[k];

                                    n_meas++;
                                }
                            }
                        }

                        // --- Energia média local ---
                        double E_local = 0.0;
                        for (int step = 0; step < t_step; ++step)
                            E_local += Energy(replicas[i]);

                        #pragma omp atomicc
                        energy_per_replica[i] += E_local / t_step;
                    }
                }
            }
        }
    }

    // --- Swaps entre réplicas ---
    for (int rep = 0; rep < rept; ++rep) {
        for (int i = 0; i < n_replicas - 1; ++i) {
            double beta_i = replicas[i].k;
            double beta_j = replicas[i + 1].k;
            double E_i = Energy(replicas[i]);
            double E_j = Energy(replicas[i + 1]);

            double delta = (beta_j - beta_i) * (E_j - E_i);
            swap_attempts++;

            if (draw_probability(gen_global) < exp(delta)) {
                std::swap(replicas[i].s, replicas[i + 1].s);
                swap_accepted++;
            }
        }
    }

    swap_acceptance_ratio = (swap_attempts > 0) ? (double)swap_accepted / swap_attempts : 0.0;

    // --- Normalização dos acumuladores ---
    for (int i = 0; i < av_s.size(); ++i)
        av_s[i] /= n_meas;

    for (int i = 0; i < av_ss.size(); ++i)
        av_ss[i] /= n_meas;

    for (int i = 0; i < n_replicas; ++i)
        energy_per_replica[i] /= rept;
}






// ===== Função paralelizada com OpenMP e sincronização dos campos =====
void parallel_tempering_multi(
    std::vector<Rede> &replicas,
    VecDoub_IO &av_s, VecDoub_IO &av_ss,
    const int t_eq, const int t_step, const int relx, const int rept,
    int n_spins, double mean, double sigma,
    const int &type, double H,
    std::vector<double> &energy_per_replica,
    double &swap_acceptance_ratio,
    std::mt19937 &gen_global
) {
    int n_replicas = replicas.size();
    int n_meas = 0;
    int swap_attempts = 0, swap_accepted = 0;

    // NÃO zerar energy_per_replica aqui! Isso é feito no loop principal se necessário.

    int target_index = 0;
    for (int i = 1; i < n_replicas; ++i)
        if (fabs(replicas[i].k - 1.0) < fabs(replicas[target_index].k - 1.0))
            target_index = i;

    int num_processors = omp_get_num_procs();
    omp_set_num_threads(num_processors);

    #pragma omp parallel
    {
        std::mt19937 gen_local(gen_global() + omp_get_thread_num());

        #pragma omp single
        {
            for (int i = 0; i < n_replicas; ++i) {
                #pragma omp task firstprivate(i)
                {
                    for (int rep = 0; rep < rept; ++rep) {
                        for (int step = 0; step < t_eq; ++step) {
                            int s_flip = draw_site(replicas[i], gen_local);
                            double dE = delta_E(replicas[i], s_flip);
                            if (dE <= 0 || draw_probability(gen_local) < exp(-replicas[i].k * dE))
                                replicas[i].s[s_flip] *= -1;
                        }

                        for (int step = 0; step < t_step; ++step) {
                            int s_flip = draw_site(replicas[i], gen_local);
                            double dE = delta_E(replicas[i], s_flip);
                            if (dE <= 0 || draw_probability(gen_local) < exp(-replicas[i].k * dE))
                                replicas[i].s[s_flip] *= -1;

                            #pragma omp critical
                            if (i == target_index && step % relx == 0) {
                                for (int j = 0; j < replicas[i].n; ++j)
                                    av_s[j] += replicas[i].s[j];
                                int ind_ss = 0;
                                for (int j = 0; j < replicas[i].n - 1; ++j)
                                    for (int k = j + 1; k < replicas[i].n; ++k)
                                        av_ss[ind_ss++] += replicas[i].s[j] * replicas[i].s[k];
                                ++n_meas;
                            }
                        }

                        double E_local = 0.0;
                        for (int step = 0; step < t_step; ++step)
                            E_local += Energy(replicas[i]);

                        #pragma omp atomic
                        energy_per_replica[i] += E_local / t_step;
                    }
                }
            }
        }
    }

    for (int rep = 0; rep < rept; ++rep) {
        for (int i = 0; i < n_replicas - 1; ++i) {
            double beta_i = replicas[i].k;
            double beta_j = replicas[i + 1].k;
            double E_i = Energy(replicas[i]);
            double E_j = Energy(replicas[i + 1]);

            double delta = (beta_j - beta_i) * (E_j - E_i);
            ++swap_attempts;

            if (draw_probability(gen_global) < exp(delta)) {
                std::swap(replicas[i].s, replicas[i + 1].s);
                ++swap_accepted;
            }
        }
    }

    swap_acceptance_ratio = (swap_attempts > 0) ? (double)swap_accepted / swap_attempts : 0.0;

    for (int i = 0; i < av_s.size(); ++i)
        av_s[i] /= n_meas;
    for (int i = 0; i < av_ss.size(); ++i)
        av_ss[i] /= n_meas;
    for (int i = 0; i < n_replicas; ++i)
        energy_per_replica[i] /= rept;
}


void wang_landau(
    Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss,
    const int t_eq, const int t_step, const int relx, const int rept,
    double f_init = 2.0, double f_min = 1e-8
) {
    unordered_map<int, double> g; // Logaritmo da densidade de estados
    unordered_map<int, int> H;    // Histograma de energias visitadas
    double f = f_init;            // Fator de modificação

    cout << "Iniciando Wang-Landau com f_init = " << f_init << " e f_min = " << f_min << endl;

    // Inicializar densidade de estados e histograma
    for (int E = -2 * r.n; E <= 2 * r.n; E += 4) {
        g[E] = 0.0;
        H[E] = 0;
    }

    // Calcular energia inicial do sistema
    int E = 0;
    for (int bond = 0; bond < r.nbonds; ++bond) {
        E += -r.J[bond] * r.s[r.no[bond * 2]] * r.s[r.no[bond * 2 + 1]];
    }
    cout << "Energia inicial: " << E << endl;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> spin_dist(0, r.n - 1);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

    while (f > f_min) {
        cout << "\nNovo ciclo com f = " << f << endl;

        for (int t = 0; t < t_step; ++t) {
            int s_flip = spin_dist(gen);
            double dE = delta_E(r, s_flip);  // Usar a mesma função de Metropolis

            cout << "Tentando flip no spin " << s_flip << " com dE = " << dE << endl;

            if (dE <= 0) {
                r.s[s_flip] = -r.s[s_flip]; // Aceita sempre se dE ≤ 0
                E += dE;
            } else {
                double p = prob_dist(gen);
                if (p < exp(g[E] - g[E + dE])) { 
                    r.s[s_flip] = -r.s[s_flip]; // Aceita com probabilidade e^(g(E) - g(E+dE))
                    E += dE;
                }
            }

            // Atualiza logaritmo da densidade de estados e histograma
            g[E] += log(f);
            H[E]++;
        }

        // Verificar achatamento do histograma
        int minH = INT_MAX, maxH = INT_MIN;
        for (const auto &entry : H) {
            if (entry.second > 0) {
                minH = min(minH, entry.second);
                maxH = max(maxH, entry.second);
            }
        }

        //cout << "Histograma: minH = " << minH << ", maxH = " << maxH << endl;

        if (minH > 0.8 * maxH) {
            f = sqrt(f);
            for (auto &entry : H) entry.second = 0; // Resetar histograma
            //cout << "Histograma achatado, reduzindo f para " << f << endl;
        }
    }

    cout << "Finalizando Wang-Landau. Iniciando amostragem...\n";

    // Processo de amostragem (igual ao Metropolis)
    for (int rep = 0; rep < rept; ++rep) {
        for (int t = 0; t < t_eq + t_step; ++t) {
            if (t >= t_eq && (int)(t - t_eq) % (int)relx == 0) {
                for (int i = 0; i < r.n; ++i) av_s[i] += r.s[i];

                int ind_ss = 0;
                for (int i = 0; i < r.n; ++i) {
                    for (int j = i; j < r.n; ++j) {
                        av_ss[ind_ss++] += r.s[i] * r.s[j];
                    }
                }
            }

            int s_flip = spin_dist(gen);
            double dE = delta_E(r, s_flip);

            if (dE <= 0) {
                r.s[s_flip] = -r.s[s_flip];
                E += dE;
            } else {
                double p = prob_dist(gen);
                if (p < exp(g[E] - g[E + dE])) {
                    r.s[s_flip] = -r.s[s_flip];
                    E += dE;
                }
            }
        }
    }

    cout << "Amostragem concluída. Normalizando valores...\n";

    // Normalizar médias
    for (int i = 0; i < av_s.size(); ++i) 
        av_s[i] /= (rept * t_step / relx);

    int ind_ss = 0;
    for (int i = 0; i < r.n; ++i) {
        for (int j = i; j < r.n; ++j) {
            av_ss[ind_ss++] /= (rept * t_step / relx);
        }
    }

    cout << "Finalização da função Wang-Landau.\n";
}


#endif
