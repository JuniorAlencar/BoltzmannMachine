from src.err_analysis import select_file, minimum_values2
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

def patterns_specific_heat(file):
    
    # Folder properties
    comp_folder = os.path.basename(os.path.dirname(file))
    # filename
    filename = os.path.basename(file)
    pattern = "a"

    if (comp_folder=="Energy"):
        pattern = r'^energia_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$'
    elif (comp_folder=="Magnetization_vs_T"):
        if (os.path.basename(file)[0:3]=="mag"):
            pattern = r'^mag_T_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$'
        else:
            pattern = r'^susc_T_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$'
    elif (comp_folder == "SpecificHeat"):
        pattern = r'^linear_specific_heat_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$'

    # Returning the name of theprops[i] originating data property, min_j, min_h, mteq, relx
    match = re.match(pattern, filename)
    variable_data = []
    if match:
        variable_data.append(match.group("filename"))
        variable_data.append(match.group("err1"))
        variable_data.append(match.group("err2"))
        variable_data.append(match.group("mteq"))
        variable_data.append(match.group("mrelx"))
        return variable_data
    else:
        print("No match found")

def load_data_specific_heat(file):
    #def load_data_specific_heat(file):
    comp_folder = os.path.dirname(os.path.dirname(file))
    
    # get ['sampleName', 'min_j_value', 'min_h_value']
    pattern = patterns_specific_heat(file)

    # filename properties
    props = ["Energy", "Magnetization", "Susceptibility", "Linear_Specific_Heat"]

    filenames = [f"{comp_folder}/Energy/energia_{pattern[0]}_err_j_{pattern[1]}_err_h_{pattern[2]}_mteq_{pattern[3]}_mrelx_{pattern[4]}.dat",
            f"{comp_folder}/Magnetization_vs_T/mag_T_{pattern[0]}_err_j_{pattern[1]}_err_h_{pattern[2]}_mteq_{pattern[3]}_mrelx_{pattern[4]}.dat",
            f"{comp_folder}/Magnetization_vs_T/susc_T_{pattern[0]}_err_j_{pattern[1]}_err_h_{pattern[2]}_mteq_{pattern[3]}_mrelx_{pattern[4]}.dat",
            f"{comp_folder}/SpecificHeat/linear_specific_heat_{pattern[0]}_err_j_{pattern[1]}_err_h_{pattern[2]}_mteq_{pattern[3]}_mrelx_{pattern[4]}.dat"]

    # create dict of dicts with propertie: {"Temperatures":Temp_values,"properties":prop_values} to all properties
    all_data = {}

    # run all properties file, saving data in all_data dictionary
    for i in range(len(filenames)):
        df = pd.read_csv(filenames[i], sep = ' ',header=None)
        df.columns = ["Temperature", props[i]]
        all_data[props[i]] = {"Temperature":df["Temperature"], props[i]:df[props[i]]}
    return props, all_data

def plotting_graph_specific_heat(Temperature, propertie, ylabel ,title):
    plt.figure(figsize=(16, 9))
    
    plt.plot(Temperature, propertie, 'o' ,color='b',ms=20, 
            alpha=0.2, mec='k')
    #plt.plot(data[prop[0]]['exp'], data[prop[0]]['exp'], color='k',linewidth=2.0)
    plt.xlabel('Temperature', fontsize=26)
    plt.ylabel(ylabel, fontsize=26)

    # Aumentar o tamanho dos ticks
    plt.tick_params(axis='both', which='major', length=10, width=2)  # Tamanho e espessura dos ticks maiores
    plt.yticks(fontsize=25)
    plt.xticks(fontsize=25)
    plt.legend(prop={"size":30}, fancybox=True, framealpha=0.0)
    plt.title(title, fontsize=30)
    
    plt.show()   
    
# Função principal para plotar os gráficos
def plotting_graphs_specific_heat(props, all_data):
    
    num_props = len(props)
    for i in range(num_props):
        # Plot each propertie
        plotting_graph_specific_heat(
        all_data[props[i]]['Temperature'], all_data[props[i]][f'{props[i]}'],props[i] ,props[i]
        )

