from src.err_analises import select_file, minimum_values
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
import pandas as pd
import re
import numpy as np
import matplotlib as mpl
from sklearn.metrics import root_mean_squared_error
mpl.rcParams['axes.linewidth'] = 1.4 #set the value globally


def pattern_properties(file):
    # desired properties
    props = ["correlation", "magnetization", "sisj", "sisjsk", "triplet"]
    # name of file
    filename = os.path.basename(file)
    # folder of propertie select
    comp_folder = os.path.basename(os.path.dirname(file))

    # file patterns for each desired property
    if(comp_folder == 'correlation'):
        pattern = r'^Pij_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)\.dat$'
    elif(comp_folder == 'magnetization'):
        pattern = r'^mag_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)\.dat$'
    elif(comp_folder == "sisj"):
        pattern = r'^sisj_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)\.dat$'
    elif(comp_folder == "sisjsk"):
        pattern = r'^sisjsk_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)\.dat$'
    elif(comp_folder == "triplet"):
        pattern = r'^Tijk_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)\.dat$'
    else:
        print(f'run properties in {props}')

    # returning the name of the originating data propertie, min_j, min_h
    match = re.match(pattern, filename)
    variable_data = []
    if match:
        variable_data.append(match.group("filename"))
        variable_data.append(match.group("err1"))
        variable_data.append(match.group("err2"))

        return variable_data

    else:
        print("No match found")

def load_data(file):
    comp_folder = os.path.dirname(os.path.dirname(file))

    # get ['sampleName', 'min_j_value', 'min_h_value']
    pattern = pattern_properties(file)

    # filename properties
    props = ["correlation", "magnetization", "sisj", "sisjsk", "triplet"]
    file_t = [comp_folder +f"/{i}/" for i in props]
    files = [f'Pij_exp_ising_{pattern[0]}_err_j_{pattern[1]}_err_h_{pattern[2]}.dat', f'mag_exp_ising_{pattern[0]}_err_j_{pattern[1]}_err_h_{pattern[2]}.dat', f'sisj_exp_ising_{pattern[0]}_err_j_{pattern[1]}_err_h_{pattern[2]}.dat',
                f'sisjsk_exp_ising_{pattern[0]}_err_j_{pattern[1]}_err_h_{pattern[2]}.dat', f'Tijk_exp_ising_{pattern[0]}_err_j_{pattern[1]}_err_h_{pattern[2]}.dat']
    filenames = [file_t[i] + files[i] for i in range(len(files))]

    # create dict of dicts with propertie: {"exp":val_exp,"ising":val_ising} to all properties
    all_data = {}

    # run all properties file, saving data in all_data dictionary
    for i in range(len(filenames)):
        df = pd.read_csv(filenames[i], sep = ' ',header=None)
        df.columns = ["exp", "ising"]
        all_data[props[i]] = {"exp":df["exp"], "ising":df["ising"]}
    
    return props, all_data

def plotting_graph(exp, ising, title):
    plt.figure(figsize=(16, 9))
    rms = root_mean_squared_error(exp, ising)
    
    plt.plot(exp, ising, 'o' ,color='b',ms=20, 
            alpha=0.2, mec='k',label=f'rmsa = {rms:.3e}')
    #plt.plot(data[prop[0]]['exp'], data[prop[0]]['exp'], color='k',linewidth=2.0)
    plt.xlabel('exp', fontsize=26)
    plt.ylabel('ising', fontsize=26)

    # Aumentar o tamanho dos ticks
    plt.tick_params(axis='both', which='major', length=10, width=2)  # Tamanho e espessura dos ticks maiores
    plt.yticks(fontsize=25)
    plt.xticks(fontsize=25)
    plt.legend(prop={"size":30}, fancybox=True, framealpha=0.0)
    plt.title(title, fontsize=30)

    plt.draw()
    xticks = plt.gca().get_xticks()
#    yticks = plt.gca().get_yticks()

    step = xticks[1]-xticks[0]
    plt.ylim([exp.min() - step, exp.max() + step])
    plt.xlim([exp.min() - step, exp.max() + step])
    straight = [i for i in exp]
    straight.insert(0,exp.min() - step)
    straight.insert(len(straight),exp.max() + step)

    plt.plot(straight, straight, color='k', linewidth=2.0)

    plt.show()

# Função principal para plotar os gráficos
def plotting_graphs(props, all_data):
    
    num_props = len(props)
    for i in range(num_props):
        # Plot each propertie
        plotting_graph(
        all_data[props[i]]['exp'], all_data[props[i]]['ising'], props[i]
        )
    