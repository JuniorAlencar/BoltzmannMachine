from src.err_analysis import select_file, minimum_values
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
import pandas as pd
import re
import numpy as np
import matplotlib as mpl
from sklearn.metrics import root_mean_squared_error
from sklearn.metrics import mean_squared_error
import glob

def pattern_properties(file):
    # Mapeamento entre pastas e prefixos dos arquivos
    property_prefixes = {
        "correlation": "Pij",
        "covariance": "Cij",
        "hi": "hi",
        "Jij": "Jij",
        "magnetization": "mag",
        "sisj": "sisj",
        "sisjsk": "sisjsk",
        "triplet": "Tijk"
    }

    # Obtém o nome do arquivo e da pasta
    filename = os.path.basename(file)
    comp_folder = os.path.basename(os.path.dirname(file))

    # Verifica se a pasta está no mapeamento
    if comp_folder not in property_prefixes:
        print(f'Unknown property folder: {comp_folder}')
        print(f'Expected one of: {list(property_prefixes.keys())}')
        return

    prefix = property_prefixes[comp_folder]

    # Cria o padrão regex dinamicamente
    pattern = (
        rf'^{prefix}_exp_ising_'
        r'(?P<filename>[a-zA-Z0-9_]+)_'
        r'err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_'
        r'err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_'
        r'mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$'
    )

    # Aplica o regex
    match = re.match(pattern, filename)
    if match:
        return [
            match.group("filename"),
            match.group("err1"),
            match.group("err2"),
            match.group("mteq"),
            match.group("mrelx")
        ]
    else:
        print(f"No match found for file: {filename}")
        return


def load_data(file):
    # Caminho até a pasta raiz das propriedades
    base_folder = os.path.dirname(os.path.dirname(file))

    # Extrai os parâmetros do nome do arquivo
    pattern = pattern_properties(file)
    if pattern is None:
        raise ValueError("Could not extract pattern from file name.")

    sample_name, err_j, err_h, mteq, mrelx = pattern

    # Mapeamento de propriedades para prefixos de arquivos e colunas
    property_map = {
        "correlation": ("Pij", "Pij"),
        "covariance": ("Cij", "Cij"),
        "hi": ("hi", "hi"),
        "Jij": ("Jij", "Jij"),
        "magnetization": ("mag", "si"),
        "sisj": ("sisj", "sisj"),
        "sisjsk": ("sisjsk", "sisjsk"),
        "triplet": ("Tijk", "Tijk"),  # Se quiser usar outro nome aqui, é só mudar
    }

    props = list(property_map.keys())
    all_data = {}

    for prop in props:
        prefix, column_prefix = property_map[prop]
        folder_path = os.path.join(base_folder, prop)
        filename = f"{prefix}_exp_ising_{sample_name}_err_j_{err_j}_err_h_{err_h}_mteq_{mteq}_mrelx_{mrelx}.dat"
        full_path = os.path.join(folder_path, filename)

        df = pd.read_csv(full_path, sep=' ')
        all_data[prop] = {
            "exp": df[f"{column_prefix}_synt"],
            "ising": df[f"{column_prefix}_ising"]
        }

    return props, all_data

def plotting_graph(exp, ising, title):
    plt.figure(figsize=(16, 9))
    rms = root_mean_squared_error(exp, ising)
    
    plt.plot(exp, ising, 'o',label=f'rmsa = {rms:.3e}')
    #plt.plot(data[prop[0]]['exp'], data[prop[0]]['exp'], color='k',linewidth=2.0)
    plt.xlabel('exp', fontsize=26)
    plt.ylabel('ising', fontsize=26)

    # Aumentar o tamanho dos ticks
    #plt.tick_params(axis='both', which='major', length=10, width=2)  # Tamanho e espessura dos ticks maiores
    #plt.yticks(fontsize=25)
    #plt.xticks(fontsize=25)
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
    

# That function clean data with ising empty and create comparative files
def clean_files():
    # desired properties
    # folder of propertie select
    folder_ising = "../Results_Metropolis/SeparateData/Cij-ising"
    
    # file patterns for each desired property
    pattern = r'^Cij_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$'
    path_in = '../Results_Metropolis/SeparateData'
    path_out = '../Results_Metropolis/Comparative'
    all_files = glob.glob(os.path.join(folder_ising, "*.dat"))
    
    for file in all_files:
        filename = os.path.basename(file)
        # Returning the name of the originating data property, min_j, min_h, mteq, relx
        match = re.match(pattern, filename)
        
        if match:
            filen = match.group("filename")
            err1 = match.group("err1")
            err2 = match.group("err2")
            mteq = match.group("mteq")
            mrelx = match.group("mrelx")
            
            # SepateData ===================>
            # ising properties
            mi_ising = path_in + f"/mi-ising/mi_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            Pij_ising = path_in + f"/Pij-ising/Pij_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            Cij_ising = path_in + f"/Cij-ising/Cij_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            sisj_ising = path_in + f"/sisj-ising/sisj_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            Tijk_ising = path_in + f"/Tijk-ising/Tijk_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            sisjsk_ising = path_in + f"/sisjsk-ising/sisjsk_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            all_files_ising = [mi_ising, Pij_ising, Cij_ising, 
                               sisj_ising, Tijk_ising,sisjsk_ising]
            # exp properties
            mi_exp = path_in + f"/mi-exp/mi_exp_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            Pij_exp = path_in + f"/Pij-exp/Pij_exp_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            Cij_exp = path_in + f"/Cij-exp/Cij_exp_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            sisj_exp = path_in + f"/sisj-exp/sisj_exp_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            Tijk_exp = path_in + f"/Tijk-exp/Tijk_exp_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            sisjsk_exp = path_in + f"/sisjsk-exp/sisjsk_exp_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            all_files_exp = [mi_exp, Pij_exp, Cij_exp, 
                               sisj_exp, Tijk_exp,sisjsk_exp]
            
            
            
            # ComparativeData ===================>
            mi_comp = f"{path_out}/magnetization/mag_exp_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            Pij_comp = f"{path_out}/correlation/Pij_exp_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            Cij_comp = f"{path_out}/covariance/Cij_exp_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            sisj_comp = f"{path_out}/sisj/sisj_exp_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            Tijk_comp = f"{path_out}/triplet/Tijk_exp_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            sisjsk_comp = f"{path_out}/sisjsk/sisjsk_exp_ising_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            all_files_comp = [mi_comp, Pij_comp, Cij_comp, 
                               sisj_comp, Tijk_comp,sisjsk_comp]
            
            all_files_rm = all_files_ising + all_files_exp + all_files_comp
            # errors
            #err = f"../Results_Metropolis/Erro/erro_{filen}_err_j_{err1}_err_h_{err2}_mteq_{mteq}_mrelx_{mrelx}.dat"
            
            # delete files empty in ising and equivalent in exp
            if os.path.exists(mi_ising) and os.path.getsize(mi_ising) == 0:
                print('deleted')
                for file_rm in all_files_rm:
                    os.remove(file_rm)
                else:
                    pass
            
            for i in range(len(all_files_ising)):
                with open(all_files_exp[i], 'r') as file_exp, open(all_files_ising[i],'r') as file_ising, open(all_files_comp[i],'w') as file_comp:
                    # Ler todas as linhas dos arquivos 1 e 2
                    lines1 = file_exp.readlines()
                    lines2 = file_ising.readlines()
                    
                    # Garantir que os arquivos tenham o mesmo número de linhas
                    if len(lines1) != len(lines2):
                        print("Os arquivos 1 e 2 não têm o mesmo número de linhas.")
                    else:
                        # Para cada linha, pegar o conteúdo dos dois arquivos e escrever no arquivo 3
                        for line1, line2 in zip(lines1, lines2):
                            # Remover espaços e quebras de linha no final de cada linha
                            col1 = line1.strip()
                            col2 = line2.strip()
                            
                            # Formatar as colunas para garantir alinhamento (exemplo: 6 casas decimais)
                            col1 = float(col1)
                            col2 = float(col2)

                            # Escrever no arquivo 3 as duas colunas, formatando com alinhamento fixo
                            file_comp.write(f"{col1:.6f} {col2:.6f}\n")

# Calculate to RMSE for each set parameters MC simulations
def RMSE_properties():
    folder_base = "../Results_Metropolis/Comparative"
    properties = ["correlation", "covariance", "magnetization", "sisj", "sisjsk","triplet"]
    patterns = [r'^Pij_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$',
            r'^Cij_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$',
            r'^mag_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$',
            r'^sisj_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$',
            r'^sisjsk_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$',
            r'^Tijk_exp_ising_(?P<filename>[a-zA-Z0-9_]+)_err_j_(?P<err1>-?\d+\.\d+e[+-]?\d+)_err_h_(?P<err2>-?\d+\.\d+e[+-]?\d+)_mteq_(?P<mteq>\d+)_mrelx_(?P<mrelx>\d+)\.dat$',]

    folders = [f"{folder_base}/{props}" for props in properties]

    folder_ising = "../Results_Metropolis/Comparative/covariance"
    dt = {"propertie":[],"sample":[],"min_j":[], "min_h":[],"teq":[],
        "relx":[],"nspins":[],"RMSE":[]}

    for j in range(len(folders)):
        all_files = glob.glob(os.path.join(folders[j], "*.dat"))

        headers = ["exp", "ising"]
        i = 0

        for file in all_files:
            filename = os.path.basename(file)
            match = re.match(patterns[j], filename)
                
            if match:
                filen = match.group("filename")
                err1 = match.group("err1")
                err2 = match.group("err2")
                mteq = match.group("mteq")
                mrelx = match.group("mrelx")
                
                if(i==0):
                    data = pd.read_csv(f"../Data/TidyData/{filen}.dat", sep=',', header=None)
                    nspins = data.shape[1]
            
                df = pd.read_csv(file, sep=" ", header=None)
                
                df.columns = headers
                
                realVals = df["exp"]
                predictedVals = df["ising"]
                mse = root_mean_squared_error(realVals, predictedVals)
                dt["propertie"].append(properties[j])
                dt["sample"].append(filen)
                dt["min_j"].append(err1)
                dt["min_h"].append(err2)
                dt["teq"].append(int(mteq)*int(nspins))
                dt["relx"].append(int(mrelx)*int(nspins))
                dt["nspins"].append(nspins)
                dt["RMSE"].append(format(mse,".2e"))
                
                i = 1

    df_final = pd.DataFrame(data=dt)
    df_final.sort_values(by="teq").reset_index(drop=True)
    return df_final

def minimum_teq(df_RMSE):
    properties = df_RMSE["propertie"].unique()
    # Criar uma lista para armazenar os DataFrames
    results_list = []
    for props in properties:
        df_f = df_RMSE[df_RMSE["propertie"]==props]
        # Selecionar os valores únicos da coluna 'teq'
        unique_teq = df_f['teq'].unique()

        for teq_value in unique_teq:
            # Filtrar o DataFrame para o valor atual de 'teq'
            subset = df_f[df_f['teq'] == teq_value]
            
            # Encontrar o índice do menor valor de 'RMSE'
            min_rmse_idx = subset['RMSE'].idxmin()
            
            # Selecionar a linha com o menor valor de 'RMSE'
            min_rmse_row = subset.loc[[min_rmse_idx]]  # Coloque entre colchetes para manter o formato de DataFrame
            
            # Adicionar a linha à lista de resultados
            results_list.append(min_rmse_row)
        
    # Concatenar todas as linhas em um único DataFrame
    results = pd.concat(results_list, ignore_index=True)
    results['order'] = results.groupby('propertie').cumcount()
    results_sorted = results.sort_values(by=['propertie', 'teq']).sort_values(by=['order']).reset_index(drop=True)
    results_sorted.drop(columns='order', inplace=True)  # Remover a coluna auxiliar 'order'
    results_sorted = results.groupby('propertie', group_keys=False).apply(lambda x: x.sort_values('teq')).reset_index(drop=True)
    results_sorted.drop(columns='order', inplace=True)
    
    return results_sorted