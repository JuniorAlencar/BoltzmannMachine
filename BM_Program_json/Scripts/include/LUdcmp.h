#ifndef LUDCMP_H
#define LUDCMP_H

struct LUdcmp
{
	Int n;
	MatDoub lu;
	VecInt indx;
	Doub d;
	LUdcmp(MatDoub_I &M);
	void solve(VecDoub_I &b, VecDoub_O &x);
	void solve(MatDoub_I &b, MatDoub_O &x);
	void inverse(MatDoub_O &ainv);
	Doub det();
	void mprove(VecDoub_I &b, VecDoub_IO &x);
	MatDoub_I &aref;
	
	void show();
};

void LUdcmp::show()
{
	int i, j;
	
	ofstream saidaLU ("LU.dat");
	
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			saidaLU << setw(11)  << lu[i][j];			
		}
		
		saidaLU << setw(11) << endl;
	}
}
	
LUdcmp::LUdcmp(MatDoub_I &M) : n(M.nrows()), lu(M), aref(M), indx(n)     //Construtor: Inicia os valores e decompoe a matriz M em LU
{
	Int i, imax, j, k;
	Doub big, temp;
	VecDoub vv(n);
	d = 1.0;
	
	for (i=0; i<n; i++)
	{
		big = 0.0;
		for (j=0; j<n; j++)
			if ((temp=abs(lu[i][j])) > big)
				big = temp;
		if (big == 0)
			throw std::runtime_error("Matriz Singular");
		vv[i] = 1.0/big;
	}
	
	for (k=0; k<n; k++)
	{
		big = 0.0;
		imax = k;
		for (i=k; i<n; i++)
		{
			temp = vv[i]*abs(lu[i][k]);
			if (temp > big)
			{
				big = temp;
				imax = i;
			}
		}
		
		if (k != imax)
		{
			for (j=0; j<n; j++)
				SWAP(lu[imax][j], lu[k][j]);
			
			d = -d;
			vv[imax] = vv[k];
		}
		
		indx[k] = imax;
		
		for (i = k +1; i<n; i++)
		{
			temp = lu[i][k] /= lu[k][k];
			for (j=k+1; j<n; j++)
				lu[i][j] -= temp*lu[k][j];
		
		}
	}
}

void LUdcmp::solve (VecDoub_I &b, VecDoub_O &x)       //Resolve o sistema de equações para um right-hand
{
	Int i, ii=0, ip, j;
	Doub sum;
	
	if (b.size() != n || x.size() != n)
		throw ("LUdcmp::solve tamanhos invalidos");
	
	for (i=0; i<n; i++)
		x[i]=b[i];
	
	for (i=0; i<n; i++)
	{
		ip = indx[i];
		sum=x[ip];
		x[ip]=x[i];
		
		if(ii != 0)
			for (j = ii - 1; j<i; j++)
				sum -= lu[i][j]*x[j];
		else
			if (sum != 0.0)
				ii = i + 1;
		x[i] = sum;
	}
	
	for ( i=n-1; i>=0; i--)
	{
		sum = x[i];
		for (j=i+1; j<n; j++)
			sum -= lu[i][j]*x[j];
		
		x[i] = sum/lu[i][i];
	}
}
						
void LUdcmp::solve (MatDoub_I &b, MatDoub_O &x)    //Resolve o sistema para varios right-hand
{
	int i, j, m=b.ncols();
	
	if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
		throw ("LUdcmp::solve tamanhos invalidos");
		
	VecDoub xx(n);
	
	for (j=0; j<m; j++)
	{
		for (i=0; i<n; i++)
			xx[i] = b[i][j];
		
		solve(xx, xx);
		
		for (i=0; i<n; i++)
			x[i][j] = xx[i];
	}
}
			
void LUdcmp::inverse(MatDoub_O &ainv)      //Calcula a matriz inversa de M
{
	int i, j;
	ainv.resize(n, n);
	
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
			ainv[i][j] = 0.0;
		
		ainv[i][i] = 1.0;
	}
	
	solve(ainv, ainv);
}

Doub LUdcmp::det()      // Determina o valor do determinante da matriz M utilizando a matriz LU
{
	int i;
	Doub dd = d;
	
	for (i=0; i<n; i++)
		dd *=lu[i][i];
		
	return dd;
}

#endif
