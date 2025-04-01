import pandas as pd
import matplotlib.pyplot as plt

# Carregar os arquivos CSV
df_uniform = pd.read_csv("uniform_vector.csv")
df_gaussian = pd.read_csv("gaussian_vector.csv")

# Criar histogramas
plt.figure(figsize=(12, 5))

# Histograma da distribuição uniforme
plt.subplot(1, 2, 1)
plt.hist(df_uniform['Uniform'], bins=20, edgecolor='black', alpha=0.7)
plt.title("Histograma - Distribuição Uniforme")
plt.xlabel("Valor")
plt.ylabel("Frequência")

# Histograma da distribuição normal (gaussiana)
plt.subplot(1, 2, 2)
plt.hist(df_gaussian['Gaussian'], bins=20, edgecolor='black', alpha=0.7)
plt.title("Histograma - Distribuição Gaussiana")
plt.xlabel("Valor")
plt.ylabel("Frequência")

# Mostrar gráficos
plt.tight_layout()
plt.show()
