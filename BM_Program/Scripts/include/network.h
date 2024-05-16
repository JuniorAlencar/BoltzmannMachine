#ifndef NETWORK_H
#define NETWORK_H



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

//Contrutor 
Rede::Rede (int m, double mmean, double ssigma, double kk, int tp, double HH) : no(2), s(m, 1.0), n(m), mean(mmean), sigma(ssigma), k(kk), type(tp), H(HH), nbonds(n*(n-1)/2), J(nbonds, 0.0), h(m, 0.0), nb(n*(n-1)), s_nb(n*(n-1)) {

	if (type == 0)
		create_bonds_random();
	else if (type == 1)
		create_bonds_tree();
		
	neighbours();
	
	double p;
	
	for (int i = 0; i < n; i++)
	{
		p = (double) rand()/RAND_MAX;
		
		if (p < 0.5)
			s[i] = -s[i];
			
	}
	
	if (H == 0)
	{
		for (int i = 0; i < n; i++)
			h[i] = 0;
	}	
	else if (H == 1)
	{
		for (int i = 0; i < n; i++)
			h[i] = (double) rand()/RAND_MAX;
	}
	else if (H == -1)
	{
		for (int i = 0; i < n; i++)
			h[i] = -(double) rand()/RAND_MAX;
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			p = (double) rand()/RAND_MAX;
		
			if (p < 0.5)
				h[i] = -(double) rand()/RAND_MAX;
			else
				h[i] = (double) rand()/RAND_MAX;
		}
	}
	
}

//Função que cria as ligaçoes aleatorias e da valores para as constante J's e h's de acordo com alguma distribuição
void Rede::create_bonds_random ()
{
	int cont = 0;
	double prob = k/ double (n-1);
	double p;
	//double JJ = gaussian();
	
	for (int i = 0; i < nbonds; i++)
	{
		p = (double) rand()/RAND_MAX;
		
		if (p <= prob)
		{
			J[i] = gaussian();
		}
	}
		
}

//Função que cria as ligações que forme uma arvore e da valores para as constantes J's e h's de acordo com alguma distribuição
void Rede::create_bonds_tree ()
{
	int cont = 0;
	double prob = k/ double (n-1);
	double p;
	VecInt root(n);
	int r0, r1;
	
	for (int i = 0; i < n; i++)
		root[i] = i;
	
	for (int i = 0; i < nbonds; i++)
	{
		p = (double) rand()/RAND_MAX;
		
		if (p <= prob)
		{
			node(i);
			
			r0 = no[0];
			r1 = no[1];
			
			while (r0 != root[r0])
				r0 = root[r0];
				
			while (r1 != root[r1])
				r1 = root[r1];
				
			if (r0 != r1)
			{
				root[r1] = r0;	
				J[i] = gaussian();
			}
		}
	}
	
	for (int i = 0; i < n; i++)
		h[i] = 0;//gaussian();
		
}


//Cria vetor de vizinhos
void Rede::neighbours()
{
	int cont = 0, lim, disl, is = 0;

	for (int i = 0; i < n; i++)
	{
		lim  = 0;
		disl = 0;
		
		for (int j = 0; j < n-1-i; j++)
		{
			nb[cont] = is;
			s_nb[cont] = j+i+1;
			lim++;
			is++;
			cont++;
		}
		
		for (int j = 0; j < n-1-lim; j++)
		{
			nb[cont] = i-1 + disl;
			s_nb[cont] = j;
			cont++;
			disl += n-2-j;
		}
	}
}

//Função que gera valores de acordo com uma ditribuição gaussiana
double Rede::gaussian ()
{

  double ymin = mean - 4.*sigma;
  double ymax = mean + 4.*sigma;
  double Pymax = 1/(sqrt (2.*M_PI)*sigma);

  // Calculate random value uniformly distributed 
  //  in range ymin to ymax
  double y = ymin + (ymax - ymin) * double (random ()) / double (RAND_MAX);

  // Calculate Py
  double Py = exp (- (y - mean) * (y - mean) / 2. / sigma / sigma) /
    sqrt (2. * M_PI) / sigma;

  // Calculate random value uniformly distributed in range 0 to Pymax
  double x = Pymax * double (random ()) / double (RAND_MAX);

  // If x > Py reject value and recalculate
  if (x > Py) 
  	return gaussian ();
  else  
  	return y;
  	
}

//Função que dado a ligação diz quais os nós
void Rede::node (int bond)
{

	int spot = 0;
	int limit = n - 1;
	int step = 1;
	int stop = 0;	
	
	while (stop == 0)
	{
		if (bond >= limit)
		{
			step++;
			limit += n - step;
			spot++;
		}
		else
			stop = 1;
	}
	
	no[0] = spot;
	no[1] = bond - (limit - n + step - 1) + spot;
	
}

#endif 
