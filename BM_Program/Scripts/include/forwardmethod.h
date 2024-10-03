#ifndef FORWARDMETHOD_H
#define FORWARDMETHOD_H

#include <random>
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

//METROPOLIS PARA BM-------------------------------------------------------------------------------
void metropolis_bm (Rede &r, VecDoub_IO &av_s, VecDoub_IO &av_ss, const int t_eq, const int t_step,  
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




#endif
