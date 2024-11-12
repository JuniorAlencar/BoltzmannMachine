import numpy as np
import pandas as pd

def load_data(file_path):
    """
    Loads and processes data from a file with comma-delimited rows into a flat numpy array.
    
    Args:
    - file_path (str): Path to the data file.
    
    Returns:
    - data_array (np.array): Flattened array of data values.
    """
    with open(file_path, 'r') as f:
        # Read each line, split by commas, and flatten into a single array
        data = []
        for line in f:
            values = line.strip().split(',')
            data.extend(map(int, values))  # Convert string values to integers
        return np.array(data)

def periodic_cross_correlation(data1, data2, max_lag=None):
    """
    Calcula a correlação cruzada entre duas amostras com condições de contorno periódicas.

    Args:
    - data1 (np.array): Array de entrada para a primeira amostra.
    - data2 (np.array): Array de entrada para a segunda amostra.
    - max_lag (int): Lag máximo para o qual calcular a correlação cruzada.

    Returns:
    - lags (np.array): Array de valores de lag.
    - cross_corr (np.array): Valores de correlação cruzada para cada lag.
    """
    n = len(data1)
    if max_lag is None:
        max_lag = n // 2  # Padrão para metade do comprimento dos dados

    lags = np.arange(-max_lag, max_lag + 1)
    cross_corr = np.zeros(len(lags))

    for i, lag in enumerate(lags):
        # Cálculo da correlação cruzada periódica usando rotação circular
        cross_corr[i] = np.mean(data1 * np.roll(data2, -lag))

    return lags, cross_corr