import os
import pandas as pd
import glob
import re

#filename: filename of data (without extension)
#j_min: minimum value Jij to MC converge (in scientific notation)
#h_min: minimum value hi to MC converge (in scientific notation)
#m_relx: multiply relx (relx = nspins*m_relx)
def multithread_local(filename, j_min, h_min, m_teq, m_relx, method):
    # file name with inputs
    file_inputs = "input_multithread.txt"
    folder = "../Scripts/"
    
    # write .txt with inputs to parallel
    # Ensure single space between values
    lst = [' '.join(map(str, (i, j, k, l, m, n))) for i, j, k, l, m, n in zip(filename, j_min, h_min, m_teq, m_relx, method)]
    
    # writing inputs
    with open(folder + file_inputs, "w") as file:
        for line in lst:
            file.write(line + "\n")
    
    # generate the outputs to run code
    
    # parallel script
    with open(folder + "parallel.sh", 'w') as a:
        # generate the outputs to run code
        a.write("#!/bin/bash\n\n")
        a.write("# Generate the executables\n")
        a.write("./compile.sh\n\n")
        a.write("# Create folders to Results if don't exist\n")
        a.write('if [ ! -d "../Results" ] || [ ! -d "../Results_Metropolis" ]; then\n')
        a.write(f"\t./CreateFolders\n")
        a.write("fi\n\n")
        a.write(f"# Run in parallel with arguments in {file_inputs}, cat read the file\n")
        a.write(f"cat {file_inputs} |parallel --colsep ' ' ./main_single.sh {{1}} {{2}} {{3}} {{4}} {{5}} {{6}}")
    
    # permission to run parallel script
    os.system("chmod +x ../Scripts/parallel.sh")
    a.close

# Load data to errors files
def load_data(file):
    df = pd.read_csv(file, sep = ' ')
    mcs_values = df['inter'].tolist()
    erroJ_values = df['erroJ'].tolist()
    erroh_values = df['erroh'].tolist()

    return mcs_values, erroJ_values, erroh_values


# find minimum values in Monte Carlo while running
def minimum_val():
    folder_erros = "../Results_Metropolis/Erro"
    pattern = r'^erro_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$'
    all_files = glob.glob(os.path.join(folder_erros,"*.dat"))
    
    data_frame = {"filename":[],"min_j_use":[],
                  "min_h_use":[],"mteq":[],"mrelx":[],"method":[]}
    
    for file in all_files:

        file_name = os.path.basename(file)
        
        mch, errJ, errH = load_data(file)
        
        min_j = format(min(errJ),".2e")
        min_h = format(min(errH),".2e")
        
        match = re.match(pattern, file_name)
        
        if match:
            data_frame["filename"].append(match.group("filename"))
            data_frame["min_j_use"].append(min_j)
            data_frame["min_h_use"].append(min_h)
            data_frame["mteq"].append(match.group("mteq"))
            data_frame["mrelx"].append(match.group("mrelx"))
            data_frame["method"].append("false")
    # DataFrame with paramers new input_multithread.txt
    df = pd.DataFrame(data=data_frame)
    df.to_csv("../Scripts/input_multithread.txt", sep=' ', header=False, index=False, mode='w+')
    