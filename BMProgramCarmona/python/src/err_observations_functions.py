import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
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
    pattern = r"err_sampleN20_j_min_([+-]?\d+\.?\d*[eE][+-]?\d+|[+-]?\d+\.?\d*)_h_min_([+-]?\d+\.?\d*[eE][+-]?\d+|[+-]?\d+\.?\d*)_t_eq_(\d+)_metropolis\.txt"

    file_name = os.path.basename(arquivo)
    match = re.match(pattern, file_name)
    
    if match:
        # Extract values as strings
        j_min_match, h_min_match, t_eq_match = match.groups()[0:3]
        
        # Convert to float and format in scientific notation
        j_min_sci = format(float(j_min_match), ".2e")
        h_min_sci = format(float(h_min_match), ".2e")
        t_eq = int(t_eq_match)  # Assuming t_eq is an integer
    
    return j_min_sci, h_min_sci, t_eq

# Função para carregar dados do arquivo TXT
def carregar_dados(arquivo, J):
    key1_values = []
    key2_values = []
    
    with open(arquivo, 'r') as f:
        for linha in f:
            valores = linha.split()  # Divide a linha pelos espaços em branco
            if len(valores) >= 2:  # Verifica se há pelo menos duas colunas
                if(J==True):
                    key1_values.append(float(valores[0]))  # Primeira coluna
                    key2_values.append(float(valores[1]))  # Segunda coluna
                else:
                    key1_values.append(float(valores[0]))  # Primeira coluna
                    key2_values.append(float(valores[2]))  # Segunda coluna
    
    return key1_values, key2_values

# Função principal para plotar o gráfico
def plotar_grafico(MCS, erro, arquivo, J):
    
    minimum = minimum_values(arquivo)
    
    if(J==True):
        ylabel = 'Err_J'
        ymin = float(minimum[0])
        label =  f' erro_min_J = {ymin}'
        label_min_err = f' erro_min_J_data = {min(erro):.2e}'
    else:
        ylabel = 'Err_h'
        ymin = float(minimum[1])
        label =  f' erro_min_h = {ymin}'
        label_min_err = f' erro_min_h_data = {min(erro):.2e}'
    
    plt.figure(figsize=(16,9))
    plt.plot(MCS, erro,color='k')
    plt.plot(MCS,ymin*np.ones(len( MCS)),'--', label=label)
    plt.xlabel('MCS',fontsize=22)
    plt.ylabel(ylabel,fontsize=22)
    plt.xlim([0, max(MCS)])
    plt.yscale("log")
    plt.title(f'Erro para t_eq = {minimum[2]}',fontsize=22)
    legend = plt.legend(prop={"size":21}, fancybox=True, framealpha=0.0)
    # Get the legend's bounding box
    legend_bbox = legend.get_window_extent()
    # Convert the legend's bounding box to figure coordinates
    fig = plt.gcf()
    legend_bbox = fig.transFigure.inverted().transform_bbox(legend_bbox)

    # Add text below the legend
    plt.text(
    x=legend_bbox.x0,  # x-coordinate relative to the figure
    y=legend_bbox.y0 - 0.05,  # y-coordinate relative to the figure (below the legend)
    s=label_min_err,  # Text to display
    transform=fig.transFigure,  # Use figure coordinates
    fontsize=22,  # Font size
    ha='left',  # Horizontal alignment
    va='top'  # Vertical alignment
    )
    plt.show()