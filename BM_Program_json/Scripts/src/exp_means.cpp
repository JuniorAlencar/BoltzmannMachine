#include "exp_means.hpp"

pair<exp_means, json> exp_calculate(string &filename, bool &use_exact){
    
    // Check if file exp exist
	string file_means = use_exact ? "../Results/" + filename + "_exp.json" : "../Results_Metropolis/" + filename + "_exp.json";
	
	// load struct with means exp values
	exp_means my_means;
	// json struct with means to append with ising means
	json json_means_exp;
	
	// If exist, just open it
	if (filesystem::exists(file_means)) {
        cout << "file exist, opening..." << endl;
		// 
		pair<exp_means, json> result = load_json_exp(file_means);
		
		return result;
    }
    else{
        // Count the number of spins in sample ---------------------------------
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
    }

    


}