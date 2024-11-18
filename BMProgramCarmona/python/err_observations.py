from src.err_observations_functions import *

# Main
if __name__ == "__main__":
    j = True
    #ss
    arquivo_selecionado = selecionar_arquivo()
    if arquivo_selecionado:
        MCS, erro = carregar_dados(arquivo_selecionado, J=j)
        plotar_grafico(MCS, erro, arquivo_selecionado ,J=j)
    else:
        print("Nenhum arquivo foi selecionado.")