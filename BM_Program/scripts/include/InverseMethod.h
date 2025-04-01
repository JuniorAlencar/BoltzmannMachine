#ifndef INVERSEDMETHOD_H
#define INVERSEMETHOD_H

#include "./LUdcmp.h"

void NaiveMeanField (Rede &r, VecDoub_I &av_s, VecDoub_I &av_ss, const double beta, const double H)
{
	int ind;

	MatDoub C(r.n, r.n, 0.0), inv_C(r.n, r.n, 0.0), J(r.n, r.n, 0.0);
	
	ind = 0;
	for (int i = 0; i < r.n; i++)
	{
		for (int j = i; j < r.n; j++)
		{
			C[i][j] = av_ss[ind] - av_s[i]*av_s[j];
			C[j][i] = C[i][j];
			ind++;
		}
	}
	
	LUdcmp nmf (C);
	
	nmf.inverse (inv_C);
	
	ind = 0;
	for (int i = 0; i < r.n-1; i++)
	{
		for (int j = i+1; j < r.n; j++)
		{
			J[i][j] = -inv_C[i][j]/beta;
			J[j][i] = J[i][j];
		
			r.J[ind] = J[i][j];
			ind++;
		}
	}
	
	if (H != 0)
	{	
		double sumJs;
	
		for (int i = 0; i < r.n; i++)
		{
			sumJs = 0;
		
			for (int j = 0; j < r.n; j++)
			{
				sumJs += J[i][j]*av_s[j];
			}
		
			r.h[i] = atanh(av_s[i])/beta - sumJs;
		}
	}
}

//-------------------------------------------------------------------------------------------------

void TAP_Equation (Rede &r, VecDoub_I &av_s, VecDoub_I &av_ss, const double beta, const double H)
{
	int ind;
	
	MatDoub C(r.n, r.n, 0.0), inv_C(r.n, r.n, 0.0), J(r.n, r.n, 0.0);
	
	ind = 0;
	for (int i = 0; i < r.n; i++)
	{
		for (int j = i; j < r.n; j++)
		{
			C[i][j] = av_ss[ind] - av_s[i]*av_s[j];
			C[j][i] = C[i][j];
			ind++;
		}
	}
	
	LUdcmp nmf (C);
	
	nmf.inverse (inv_C);
	
	if (H == 0)
	{
		ind = 0;
		for (int i = 0; i < r.n-1; i++)
		{
			for (int j = i+1; j < r.n; j++)
			{
				J[i][j] = -inv_C[i][j]/beta;
				J[j][i] = J[i][j];
		
				r.J[ind] = J[i][j];
				ind++;
			}
		}
	}
	else
	{	
		ind = 0;
		for (int i = 0; i < r.n-1; i++)
		{
			for (int j = i+1; j < r.n; j++)
			{
				J[i][j] = (sqrt(1-8*av_s[i]*av_s[j]*inv_C[i][j]) - 1)/(4*beta*av_s[i]*av_s[j]);
				J[j][i] = J[i][j];
		
				r.J[ind] = J[i][j];
				ind++;
			}
		}
	
		double sumJs, sumJ2s;
	
		for (int i = 0; i < r.n; i++)
		{
			sumJs = 0;
		
			for (int j = 0; j < r.n; j++)
			{
				sumJs += J[i][j]*av_s[j];
				sumJ2s += pow(J[i][j], 2)*(1-pow(av_s[j], 2));
			}
		
			r.h[i] = atanh(av_s[i])/beta - sumJs + beta*av_s[i]*sumJ2s;
		}
	}

}

#endif
