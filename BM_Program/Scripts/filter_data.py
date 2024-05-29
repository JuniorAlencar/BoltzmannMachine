
import numpy as np
import pandas as pd

#n_rows: number of samples
#n_colums: number of spins
#filename: number of file with data bin

#dataTeste
def filter_data():
    # Get user input
    filename = input("Enter the filename (without extension): ")
    n_rows = int(input("Enter the number of rows to select(num_samples): "))
    n_colums = int(input("Enter the number of columns to select(num_spins): "))
    
    df = pd.read_csv(f"../Data/TidyData/{filename}.dat",sep=',', header=None)

    column_names= [f's_{i}' for i in range(df.shape[1])]
    df.columns = column_names
    df_filter = df.iloc[:n_rows, :n_colums]
    df_filter.to_csv(f"../Data/TidyData/{filename}_filter.dat",index=False, header=None)

filter_data()