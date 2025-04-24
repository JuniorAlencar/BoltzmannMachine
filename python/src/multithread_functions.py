import os
import pandas as pd
import glob
import re

#filename: filename of data (without extension)
#j_min: minimum value Jij to MC converge (in scientific notation)
#h_min: minimum value hi to MC converge (in scientific notation)
#m_relx: multiply relx (relx = nspins*m_relx)
def multithread_local(filename, j_min, h_min, m_teq, m_relx, seed, method):
    # file name with inputs
    file_inputs = "input_multithread.txt"
    folder = "../shells/"
    
    # write .txt with inputs to parallel
    # Ensure single space between values
    lst = [' '.join(map(str, (i, j, k, l, m, n, o))) for i, j, k, l, m, n, o in zip(filename, j_min, h_min, m_teq, m_relx, seed, method)]
    
    # writing inputs
    with open(folder + file_inputs, "w") as file:
        for line in lst:
            file.write(line + "\n")
    
    # generate the outputs to run code
    
    # parallel script
    with open(folder + "parallel.sh", 'w') as a:
        # generate the outputs to run code
        a.write("#!/bin/bash\n\n# Generate the executables\n./compile.sh\n\n# Create folders to Results if don't exist\n")
        a.write('if [ ! -d "../Results" ] || [ ! -d "../Results_Metropolis" ]; then\n')
        a.write(f"\t./CreateFolders\n")
        a.write("fi\n\n")
        a.write(f"# Run in parallel with arguments in {file_inputs}, cat read the file\n")
        a.write(f"cat {file_inputs} |parallel --colsep ' ' ./BoltmannMachine.sh {{1}} {{2}} {{3}} {{4}} {{5}} {{6}} {{7}}")
    
    # permission to run parallel script
    os.system("chmod +x ../shells/parallel.sh")
    a.close

# Load data to errors files
def load_data(file):
    df = pd.read_csv(file, sep = ' ')
    mcs_values = df['inter'].tolist()
    erroJ_values = df['erroJ'].tolist()
    erroh_values = df['erroh'].tolist()

    return mcs_values, erroJ_values, erroh_values


# find minimum values in Monte Carlo while running
def minimum_val(filename, method, comb_parms):
    """
        filename (string): name of file (without extension .dat)
        method (string): name of algorithm used (metropolis, exact ou parallel_tempering)
        comb_parms (tuple): all combinations of (m_relx, m_teq)
    
    """
    folder_erros = f"../Results/{method}/Erro"
    pattern = fr'^erro_{re.escape(filename)}_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)_seed_(?P<seed>\d+)\.dat$'
    all_files = glob.glob(os.path.join(folder_erros, "*.dat"))

    data_frame = {
        "filename": [],
        "min_j_use": [],
        "min_h_use": [],
        "mteq": [],
        "mrelx": [],
        "seed": [],
        "method": []
    }

    for file in all_files:
        file_name = os.path.basename(file)
        match = re.match(pattern, file_name)

        if match:
            # Pega os valores e converte para int
            mteq_val = int(match.group("mteq"))
            mrelx_val = int(match.group("mrelx"))
            seed_val  = int(match.group("seed"))
            min_j_used = float(match.group("err1"))
            min_h_used = float(match.group("err1"))

            # Verifica se o par está em comb_parms
            if (mteq_val, mrelx_val) in comb_parms:
                mch, errJ, errH = load_data(file)
                
                # Minimum from file
                min_j = format(min(errJ), ".2e")
                min_h = format(min(errH), ".2e")
                
                # Minimum in simulation

                # If system converged, ignore file
                if (float(min_j) <= min_j_used and float(min_h) <= min_h_used): 
                    pass
                else:
                    data_frame["filename"].append(filename)
                    data_frame["min_j_use"].append(min_j)
                    data_frame["min_h_use"].append(min_h)
                    data_frame["mteq"].append(mteq_val)
                    data_frame["mrelx"].append(mrelx_val)
                    data_frame["method"].append(method)
                    data_frame["seed"].append(seed_val)
    # DataFrame with paramers new input_multithread.txt
    df = pd.DataFrame(data=data_frame)
    df.to_csv("../shells/input_multithread.txt", sep=' ', header=False, index=False, mode='w+')
    