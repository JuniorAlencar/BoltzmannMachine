import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
import pandas as pd
import re
import numpy as np
import matplotlib as mpl
from matplotlib.animation import FuncAnimation
import time

mpl.rcParams['axes.linewidth'] = 1.4 #set the value globally

# Função para selecionar o file .dat
def select_file():
    root = tk.Tk()
    root.withdraw()  # Esconde a janela principal do Tkinter
    
    
    file = filedialog.askopenfilename(
        title="Selecione um arquivo TXT",
        filetypes=(("Text files", "*.dat"), ("All files", "*.*"))
    )
    
    return file

def minimum_values(file):
    file_name = os.path.basename(file)
    
    if(file_name[5:14] == "sampleN30"):
        pattern = r"erro_sampleN30_err_j_([-+]?[0-9]*\.?[0-9]+[eE][-+]?[0-9]+)_err_h_([-+]?[0-9]*\.?[0-9]+[eE][-+]?[0-9]+)\.dat"
    elif(file_name[5:14] == "sampleN20"):
        pattern = r"erro_sampleN20_err_j_([-+]?[0-9]*\.?[0-9]+[eE][-+]?[0-9]+)_err_h_([-+]?[0-9]*\.?[0-9]+[eE][-+]?[0-9]+)\.dat"
    else:
        pattern = r"erro_sampleCarmona_err_j_([-+]?[0-9]*\.?[0-9]+[eE][-+]?[0-9]+)_err_h_([-+]?[0-9]*\.?[0-9]+[eE][-+]?[0-9]+)\.dat"
    match = re.match(pattern, file_name)
    
    if match:
        # Extract values as strings
        j_min_match, h_min_match = match.groups()
        
        # Convert to float and format in scientific notation
        j_min_sci = format(float(j_min_match), ".2e")
        h_min_sci = format(float(h_min_match), ".2e")
        
    return j_min_sci, h_min_sci

def minimum_values2(file):
    file_name = os.path.basename(file)
    
    if(file_name[5:14] == "sampleN30"):
        pattern = r"erro_sampleN30_err_j_([-+]?\d+\.\d+[eE][-+]?\d+)_err_h_([-+]?\d+\.\d+[eE][-+]?\d+)_mteq_(\d+)_mrelx_(\d+)\.dat"
    elif(file_name[5:14] == "sampleN20"):
        pattern = r"erro_sampleN20_err_j_([-+]?\d+\.\d+[eE][-+]?\d+)_err_h_([-+]?\d+\.\d+[eE][-+]?\d+)_mteq_(\d+)_mrelx_(\d+)\.dat"
    else:
        pattern = r"erro_sampleCarmona_err_j_([-+]?\d+\.\d+[eE][-+]?\d+)_err_h_([-+]?\d+\.\d+[eE][-+]?\d+)_mteq_(\d+)_mrelx_(\d+)\.dat"
    match = re.match(pattern, file_name)
    
    if match:
        # Extract values as strings
        j_min_match, h_min_match, t_eq, relx = match.groups()
        
        # Convert to float and format in scientific notation
        #j_min_sci = format(float(j_min_match), ".2e")
        #h_min_sci = format(float(h_min_match), ".2e")
        
    return j_min_match, h_min_match, t_eq, relx


# Função para carregar dados do file DAT
def load_data(file):
    df = pd.read_csv(file, sep = ' ')
    mcs_values = df['inter'].tolist()
    erroJ_values = df['erroJ'].tolist()
    erroh_values = df['erroh'].tolist()

    return mcs_values, erroJ_values, erroh_values

def plotting_graph(mcs, erro, ylabel, ymin, label, label_min_err):
    plt.figure(figsize=(16, 9))
    plt.plot(mcs, erro, color='k')
    plt.plot(mcs, ymin * np.ones(len(mcs)), '--', label=label, linewidth=1.4)
    plt.xlabel('MCS', fontsize=22)
    plt.ylabel(ylabel, fontsize=22)
    plt.xlim([0, max(mcs)])
    plt.yscale("log")
    plt.title('Erro para t_eq = 3000', fontsize=22)
    
    legend = plt.legend(prop={"size": 21}, fancybox=True, framealpha=0.0)
    legend_bbox = legend.get_window_extent()
    
    fig = plt.gcf()
    legend_bbox = fig.transFigure.inverted().transform_bbox(legend_bbox)
    
    # Ajustando o posicionamento do texto para ficar abaixo da legenda
    plt.text(
        x=legend_bbox.x0, 
        y=legend_bbox.y0 - 0.05, 
        s=label_min_err, 
        transform=fig.transFigure, 
        fontsize=22, 
        ha='left', 
        va='top'
    )
    
    plt.show()

def plotting_graph2(mcs, erro, ylabel, ymin, label, label_min_err, t_eq, relx):
    plt.figure(figsize=(16, 9))
    plt.plot(mcs, erro, color='k')
    plt.plot(mcs, ymin * np.ones(len(mcs)), '--', label=label, linewidth=1.4)
    plt.xlabel('MCS', fontsize=22)
    plt.ylabel(ylabel, fontsize=22)
    plt.xlim([0, max(mcs)])
    plt.yscale("log")
    plt.title(f'Erro para t_eq = {t_eq}, relx = {relx}', fontsize=22)
    
    legend = plt.legend(prop={"size": 21}, fancybox=True, framealpha=0.0)
    legend_bbox = legend.get_window_extent()
    
    fig = plt.gcf()
    legend_bbox = fig.transFigure.inverted().transform_bbox(legend_bbox)
    
    # Ajustando o posicionamento do texto para ficar abaixo da legenda
    plt.text(
        x=legend_bbox.x0, 
        y=legend_bbox.y0 - 0.05, 
        s=label_min_err, 
        transform=fig.transFigure, 
        fontsize=22, 
        ha='left', 
        va='top'
    )
    
    plt.show()
#mcs, erro, ylabel, ymin, label, label_min_err, t_eq, relx
# Função principal para plotar os gráficos
def plotting_graphs2(mcs, erroJ, erroh, file):
    minimum = minimum_values2(file)
    
    # Plot para Err_J
    plotting_graph2(
        mcs, erroJ, 'Err_J', float(minimum[0]), 
        f'erro_min_J = {float(minimum[0]):.2e}', 
        f'erro_min_J_data = {min(erroJ):.2e}',minimum[2], minimum[3]
    )
    
    # Plot para Err_h
    plotting_graph2(
        mcs, erroh, 'Err_h', float(minimum[1]), 
        f'erro_min_h = {float(minimum[1]):.2e}', 
        f'erro_min_h_data = {min(erroh):.2e}',minimum[2], minimum[3]
    )
    print(f'data_min_h = {min(erroh):.2e},data_min_j = {min(erroJ):.2e}, t_eq = {minimum[2]}, relx =  {minimum[3]}')

# Função principal para plotar os gráficos
def plotting_graphs(mcs, erroJ, erroh, file):
    minimum = minimum_values(file)
    
    # Plot para Err_J
    plotting_graph(
        mcs, erroJ, 'Err_J', float(minimum[0]), 
        f'erro_min_J = {float(minimum[0]):.2e}', 
        f'erro_min_J_data = {min(erroJ):.2e}'
    )
    
    # Plot para Err_h
    plotting_graph(
        mcs, erroh, 'Err_h', float(minimum[1]), 
        f'erro_min_h = {float(minimum[1]):.2e}', 
        f'erro_min_h_data = {min(erroh):.2e}'
    )
    print(f'data_min_h = {min(erroh):.2e},data_min_j = {min(erroJ):.2e}')


# Function to load the data
def load_data_interactive(file):
    try:
        df = pd.read_csv(file, sep=r'\s+', dtype={'inter': int, 'erroJ': float, 'erroh': float})
        return df['inter'], df['erroJ'], df['erroh']
    except Exception as e:
        print(f"Error reading file: {e}")
        return [], [], []
    
# Function to update the plots
def update_plot(file_path):
    # Load the data
    inter, erroJ, erroh = load_data_interactive(file_path)
    
    # Initialize the first plot (ErroJ)
    fig1, ax1 = plt.subplots(figsize=(16, 9))
    ax1.set_title('ErroJ over Time', fontsize=22)
    ax1.set_xlabel('MCS', fontsize=18)
    ax1.set_ylabel('ErroJ (log scale)', fontsize=18)
    ax1.set_yscale('log')
    line1, = ax1.plot([], [], color='blue', label='ErroJ')
    ax1.legend(fontsize=14)

    # Initialize the second plot (Erroh)
    fig2, ax2 = plt.subplots(figsize=(16, 9))
    ax2.set_title('Erroh over Time', fontsize=22)
    ax2.set_xlabel('MCS', fontsize=18)
    ax2.set_ylabel('Erroh (log scale)', fontsize=18)
    ax2.set_yscale('log')

    line2, = ax2.plot([], [], color='red', label='Erroh')
    ax2.legend(fontsize=14)

    # Real-time plotting loop
    print("Press Ctrl+C to stop the real-time plot.")
    try:
        while True:
            if os.path.exists(file_path):
                inter, erroJ, erroh = load_data_interactive(file_path)
                if len(inter) > 0:
                    # Update ErroJ plot
                    line1.set_data(inter, erroJ)
                    ax1.set_xlim(0, max(inter))
                    ax1.set_ylim(min(erroJ) * 0.9, max(erroJ) * 1.1)
                    fig1.canvas.draw_idle()

                    # Update Erroh plot
                    line2.set_data(inter, erroh)
                    ax2.set_xlim(0, max(inter))
                    ax2.set_ylim(min(erroh) * 0.9, max(erroh) * 1.1)
                    fig2.canvas.draw_idle()

                plt.pause(1)  # Pause for a moment to refresh the plots
            time.sleep(1)  # Check for updates every second
    except KeyboardInterrupt:
        print("Real-time plot stopped.")
        plt.close(fig1)
        plt.close(fig2)