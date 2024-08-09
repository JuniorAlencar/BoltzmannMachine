import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import sys
#from matplotlib.animation import FuncAnimation
import os

# Abre o arquivo cpp do BMFinal e retorna o min_erro_j e min_erro_h
def extract_formatted_values_from_cpp():
    file_path = "../Scripts/BMfinal.cpp"
    results = {
        "min_erro_j": None,
        "min_erro_h": None
    }

    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

        for line in lines:
            # Search for pattern double min_erro_j
            match_j = re.search(r'double\s+min_erro_j\s*=\s*(-?\d+\.\d+e[-+]?\d+)', line)
            if match_j:
                value_j = float(match_j.group(1))
                results["min_erro_j"] = format(value_j, ".1e")  # Format as scientific notation with minimal precision

            # Search for pattern double min_erro_h
            match_h = re.search(r'double\s+min_erro_h\s*=\s*(-?\d+\.\d+e[-+]?\d+)', line)
            if match_h:
                value_h = float(match_h.group(1))
                results["min_erro_h"] = format(value_h, ".1e")  # Format as scientific notation with minimal precision

            # Break the loop if both patterns are found
            if results["min_erro_j"] is not None and results["min_erro_h"] is not None:
                break

    return results

def erro_parms(parms,N_spins, save):
    # Create folder to test_erro
    folder = f"./tests_erro/{parms}"
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    min_values = extract_formatted_values_from_cpp()

    df = pd.read_csv(f"../Results/Erro/erro_sampleN{N_spins}.dat",sep=' ',header=None)
    df.columns = ["MCS", "Erro_J", "Erro_h"]
    x = df["MCS"]
    plt.figure(figsize=(16,9))
    if parms=="J":
        y = df["Erro_J"]
        if len(y)==1:
            return print(f'erro_{parms} Atinge mínimo com um MCH para Nspins = {N_spins}')
        else:
            y_lim = np.ones(len(y))*9.712739e-04
    
    elif parms=="h":
        y = df["Erro_h"]
        if len(y)==1:
            return print(f'erro_{parms} Atinge mínimo com um MCH para Nspins = {N_spins}')
        else:
            y_lim = np.ones(len(y))*6.670033e-04
    
    # Clear previous plot
    #print(float(min_values['min_erro_j']))
    plt.yscale("log")
    plt.plot(x,y,'-',color='k')
    plt.plot(x,y_lim,'--',color='red')
    plt.xlabel("MCS",fontsize=22)
    plt.ylabel(f"Erro {parms}",fontsize=22)
    plt.xlim(min(x),max(x))
    plt.title(f"Erro {parms} para Nspins = {N_spins}",fontsize=25)
    
    if(save == True):
        j_min = min_values["min_erro_j"]
        h_min = min_values["min_erro_h"]
        
        plt.savefig(folder + f"/sampleN{N_spins}_parm_{parms}_j_{j_min}_h_{h_min}.pdf", dpi=300)
        
    plt.show()

def erro_min(N_spins, value):
    df = pd.read_csv(f"../Results/Erro/erro_sampleN{N_spins}.dat",sep=' ',header=None)
    df.columns = ["MCS", "Erro_J", "Erro_h"]

    err_min_J = df["Erro_J"].min()
    err_min_idx_J = df["Erro_J"].idxmin()
    err_min_j="{:e}".format(err_min_J)

    err_min_H = df["Erro_h"].min()
    err_min_idx_h = df["Erro_h"].idxmin()
    err_min_h="{:e}".format(err_min_H)
    if(value==True):
        return print(err_min_j, err_min_h)
    else:
        return print(df["Erro_h"][76363])
        #return print(err_min_idx_J, err_min_idx_h)
    # if parms=="J":
    #     err_min_J = df["Erro_J"].min()
    #     err_min_idx_J = df["Erro_J"].idxmin()
    #     err_min_j="{:e}".format(err_min_J)
        
    #     return print(err_min_j, err_min_idx_J)
    
    # elif parms=="h":
    #     err_min_H = df["Erro_h"].min()
    #     err_min_idx_h = df["Erro_h"].idxmin()
    #     err_min_h="{:e}".format(err_min_H)
        
    #     return print(err_min_h, err_min_idx_h)

# nspins: number of spins, parms: h or J
if __name__ == "__main__":
    nspins = int(sys.argv[1])
    parms = str(sys.argv[2])
    save = sys.argv[3].lower() == 'true'
    #erro_min(nspins,True)
    erro_parms(parms, nspins, save)