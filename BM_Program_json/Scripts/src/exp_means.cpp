#include "exp_means.hpp"

exp_means exp_mean_calculate::exp_calculate(const string &filename){
    // Add flag class to use all functions defined
    js_funct rnd;
    
    // Check if file exp exist
	string file_means = "../Results/" + filename + "_exp.json";
	
	// load struct with means exp values
	exp_means my_means;
	
	// If exist, just open it
	if (filesystem::exists(file_means)) {
        cout << "file exist, opening..." << endl;
		// 
		my_means = rnd.load_json_exp(filename);
    }
    else{
        // Count the number of spins in sample ---------------------------------
        cout << "file experimental don't exist, calculate..." << endl;
        
        string file_input = "../Data/TidyData/" + filename + ".dat";
        ifstream data_input (file_input.c_str());
        string first_line;
        string a;

        getline(data_input, first_line);
        int N = count(first_line.begin(), first_line.end(), ',');
        
        // Number of spins
        int nspins = N + 1; 
        
        // count the number of samples
        int m = 0;
        
        int number = N + 1;
        
        while (getline(data_input, first_line))
        {
            m++;
        }

        m++;

        cout << "Number of spins = " << number << endl;
        cout << "Number of samples = " << m << endl;

        // return to the beginning of the file
        data_input.clear();
        data_input.seekg (0, ios::beg);

        
        // Create the matrix to receive the values
        vector<vector<int>> M(m, vector<int>(N));
        
        //getline(data_input, first_line);

        for (int w = 0; w < m; w++)
        {
            for(int p = 0; p < N; p++)
            {
                if (p < N - 1)
                {
                    getline(data_input, a, ',');
                }
                else
                {
                    getline(data_input, a);
                }           
                // Update matrix
                M[w][p] = stoi(a);				
            }

        }

        // Close file
        data_input.close();
        
        // Number of combinations in product (xi*xj*xk)
        int n_triplet = N*(N-1)*(N-2)/6;
	    // Number of combinations in product (xi*xj)
        int n_duplet  = N*(N-1)/2;
        
		// First, second and third moment
		vector<double> s(N, 0.0), ss(n_duplet, 0.0), sss(n_triplet, 0.0);  
		// Auxiliary matrix to calculate the triplet and third moment
		vector<vector<double>> M_ss(N, vector<double>(N));
		
		// First moment (magnetization)
		for (int p = 0; p < N; p++)
		{
			for(int w = 0; w < m; w++)
			{
				s[p] += M[w][p];
			}
			
			s[p] /= m;
		}

		// Second moment
		int ind = 0;
		for (int p = 0; p < N-1; p++)
		{
			for(int pp = p+1; pp < N; pp++)
			{		
				for (int w = 0; w < m; w++)
				{
					ss[ind] += M[w][p]*M[w][pp];
				}
			
				ss[ind] /= m;

				M_ss[p][pp] = ss[ind];
				M_ss[pp][p] = M_ss[p][pp];
			
				ind++;
			}
		}
		// Third moment
		ind = 0;
		for (int i = 0; i < N-2; i++)
		{
			for (int j = i+1; j < N-1; j++)
			{
				for (int k = j+1; k < N; k++)
				{
					for (int s = 0; s < m; s++)
					{
						sss[ind] += M[s][i]*M[s][j]*M[s][k];
					}
					
					sss[ind] /= m;
					ind++;
				}
			}
		}

		// Triplet
		vector<double> Triplet(n_triplet, 0.0);
		
		ind = 0;
		for (int i = 0; i < N-2; i++)
		{
			for (int j = i+1; j < N-1; j++)
			{
				for (int k = j+1; k < N; k++)
				{
					Triplet[ind] = sss[ind] - s[i]*M_ss[j][k] - s[j]*M_ss[i][k] 
									- s[k]*M_ss[i][j] + 2*s[i]*s[j]*s[k];
					
									
					ind++;
				}
			}
		}

		// Covariance experimental (Cij)
		vector<double> C(N*(N-1)/2, 0.0);
		
		// Correlation experimental (Pij)
		vector<double> Pearson(n_duplet);

			ind = 0;
		for (int p = 0; p < N-1; p++)
		{
			for (int pp = p+1; pp < N; pp++)
			{
				C[ind] = ss[ind] - s[p]*s[pp];

				Pearson[ind] = C[ind]/sqrt((1 - pow(s[p], 2))*(1 - pow(s[pp], 2)));
				
				ind++;
			}
		}
		
		// Update struct exp_means with values
		my_means.av_s = s;
		my_means.av_ss = ss;
		my_means.av_sss = sss;
		my_means.Cij_exp = C;
		my_means.Pij_exp = Pearson;
		my_means.Tijk_exp = Triplet;
		
		// Json saved and created
		rnd.create_json_exp(my_means, filename);
    }
    
    return my_means;
}