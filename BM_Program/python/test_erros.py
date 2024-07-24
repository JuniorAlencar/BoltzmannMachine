import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import sys
from matplotlib.animation import FuncAnimation

# Abre o arquivo cpp do BMFinal e retorna o min_erro_j e min_erro_h
# def extract_formatted_values_from_cpp():
#     file_path = "../Scripts/BMfinal.cpp"
#     results = {
#         "min_erro_j": None,
#         "min_erro_h": None
#     }

#     with open(file_path, 'r', encoding='utf-8') as f:
#         lines = f.readlines()

#         for line in lines:
#             # Search for pattern double min_erro_j
#             match_j = re.search(r'double\s+min_erro_j\s*=\s*(-?\d+\.\d+e[-+]?\d+)', line)
#             if match_j:
#                 value_j = float(match_j.group(1))
#                 results["min_erro_j"] = format(value_j, ".1e")  # Format as scientific notation with minimal precision

#             # Search for pattern double min_erro_h
#             match_h = re.search(r'double\s+min_erro_h\s*=\s*(-?\d+\.\d+e[-+]?\d+)', line)
#             if match_h:
#                 value_h = float(match_h.group(1))
#                 results["min_erro_h"] = format(value_h, ".1e")  # Format as scientific notation with minimal precision

#             # Break the loop if both patterns are found
#             if results["min_erro_j"] is not None and results["min_erro_h"] is not None:
#                 break

#     return results



def erro_parms(parms,N_spins):
    #min_values = extract_formatted_values_from_cpp()

    df = pd.read_csv(f"../Results/Erro/erro_sampleN{N_spins}.dat",sep=' ',header=None)
    df.columns = ["MCS", "Erro_J", "Erro_h"]
    x = df["MCS"]

    if parms=="J":
        y = df["Erro_J"]
        if len(y)==1:
            return print(f'erro_{parms} Atinge mínimo com um MCH para Nspins = {N_spins}')
        else:
            y_lim = np.ones(len(y))*df["Erro_J"].min()
            #y_lim = np.ones(len(y))*float(min_values['min_erro_j'])
    
    elif parms=="h":
        y = df["Erro_h"]
        if len(y)==1:
            return print(f'erro_{parms} Atinge mínimo com um MCH para Nspins = {N_spins}')
        else:
            y_lim = np.ones(len(y))*df["Erro_h"].min()
    
    # Clear previous plot
    #print(float(min_values['min_erro_j']))
    plt.yscale("log")
    plt.plot(x,y,'-',color='k')
    plt.plot(x,y_lim,'--',color='red')
    plt.xlabel("MCS",fontsize=22)
    plt.ylabel(f"Erro {parms}",fontsize=22)
    plt.xlim(min(x),max(x))
    plt.title(f"Erro {parms} para Nspins = {N_spins}",fontsize=25)
    plt.show()

# nspins: number of spins, variable: propertie
if __name__ == "__main__":
    nspins = int(sys.argv[1])
    parms = str(sys.argv[2])
    erro_parms(parms,nspins)
    #test_parms(variable,nspins)    
# Initialize the plot
#fig, ax = plt.subplots()

# Create an animation
#ani = FuncAnimation(fig, erro_parms(parms,N_spins), frames=np.arange(0, 10), interval=1000)  # Update every 10 seconds (10000 milliseconds)

# Show the plot
#plt.show()