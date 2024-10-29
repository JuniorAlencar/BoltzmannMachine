import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
import pandas as pd
import re
import numpy as np

# Função para selecionar o arquivo .txt
def selecionar_arquivo():
    root = tk.Tk()
    root.withdraw()  # Esconde a janela principal do Tkinter
    
    
    arquivo = filedialog.askopenfilename(
        title="Selecione um arquivo TXT",
        filetypes=(("Text files", "*.txt"), ("All files", "*.*"))
    )
    
    return arquivo

def minimum_values(arquivo):
    # Sample filename
    file_name = os.path.basename(arquivo)

    # Pattern to err_sample if N=20 and N=30    
    if(file_name[5:14] == "sampleN30"):
        pattern = r"err_sampleN30_j_min_([+-]?\d+\.?\d*[eE][+-]?\d+|[+-]?\d+\.?\d*)_h_min_([+-]?\d+\.?\d*[eE][+-]?\d+|[+-]?\d+\.?\d*)_t_eq_(\d+)\.txt"
    else:
        pattern = r"err_sampleN20_j_min_([+-]?\d+\.?\d*[eE][+-]?\d+|[+-]?\d+\.?\d*)_h_min_([+-]?\d+\.?\d*[eE][+-]?\d+|[+-]?\d+\.?\d*)_t_eq_(\d+)\.txt"

    # Match the pattern
    match = re.match(pattern, file_name)

    if match:
        # Extract values as strings
        value1, value2, value3 = match.groups()
        
        try:
            # Convert to float and format in scientific notation
            j_min = format(float(value1), ".2e")
            h_min = format(float(value2), ".2e")
            t_eq = format(float(value3), ".2e")  # Assuming value3 is an integer
            
            return j_min, h_min, t_eq
        except ValueError as e:
            print(f"Error converting value: {e}")
    else:
        print("No match found")

# Função para carregar dados do arquivo TXT
def carregar_dados(arquivo):
    df = pd.read_csv(arquivo, sep = ' ')
    
    mcs_values = df['iter'].tolist()
    erroJ_values = df['err_J'].tolist()
    erroh_values = df['err_h'].tolist()

    return mcs_values, erroJ_values, erroh_values

def plotar_grafico(mcs, erro, ylabel, ymin, label, label_min_err, title):
    plt.figure(figsize=(16, 9))
    plt.plot(mcs, erro, color='k')
    plt.plot(mcs, ymin * np.ones(len(mcs)), '--', label=label, linewidth=1.4)
    plt.xlabel('MCS', fontsize=22)
    plt.ylabel(ylabel, fontsize=22)
    plt.title(title)
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

# Função principal para plotar os gráficos
def plotar_graficos(mcs, erroJ, erroh, arquivo):
    minimum = minimum_values(arquivo)
    
    # Plot para Err_J
    plotar_grafico(
        mcs, erroJ, 'Err_J', float(minimum[0]), 
        f'erro_min_J = {float(minimum[0]):.2e}', 
        f'erro_min_J_data = {min(erroJ):.2e}',
        f't_eq = {minimum[2]}'
    )
    
    # Plot para Err_h
    plotar_grafico(
        mcs, erroh, 'Err_h', float(minimum[1]), 
        f'erro_min_h = {float(minimum[1]):.2e}', 
        f'erro_min_h_data = {min(erroh):.2e}',
        f't_eq = {minimum[2]}'
    )