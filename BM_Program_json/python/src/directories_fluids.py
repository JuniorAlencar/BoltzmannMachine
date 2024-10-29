import os
import glob

# Função para encontrar arquivos .dat dentro de pastas "configs"
def find_dat_files_in_configs(root_dir, output_file):
    dat_files = []
    
    if os.path.isfile(output_file):
        print("file_exist")
        pass
    
    else:
        # Percorrer toda a estrutura de diretórios
        for root, dirs, files in os.walk(root_dir):
            # Verificar se estamos numa pasta "configs"
            if 'configs' in root.lower():
                for file in files:
                    # Checar se o arquivo termina com .dat
                    if file.endswith('.dat'):
                        # Armazenar o path completo do arquivo
                        dat_files.append(os.path.join(root, file))

        # Escrever o path dos arquivos .dat em um arquivo de texto
        with open(output_file, 'w') as f:
            for file_path in dat_files:
                f.write(file_path + '\n')

        print(f"Arquivo de saída salvo em: {output_file}")
    
    
def copy_directory_structure(src_dir, dest_dir):
    if os.path.isdir("../Results/Flow_Spin"):
        print("structure folders exist")
        pass
    else:
        # Iterar sobre todos os diretórios e subdiretórios
        for root, dirs, files in os.walk(src_dir):
            # Criar o path correspondente no diretório de destino
            relative_path = os.path.relpath(root, src_dir)  # path relativo da estrutura
            dest_path = os.path.join(dest_dir, relative_path)  # path de destino

            # Criar o diretório no destino
            os.makedirs(dest_path, exist_ok=True)
            print(f"Diretório criado: {dest_path}")

# Função para pegar o path após a segunda barra
def second_bar(path):
    partes = path.split('/')
    return '/' + '/'.join(partes[3:]) if len(partes) > 3 else path

def binning_data():
    # create executable(if dont exist) to binning data in c++ -------------------------
    exec = "./bin_data"
    if os.path.isfile(exec):
        pass
    else:
        os.system("g++ -std=c++17 -o bin_data ../src/clean_data.cpp")
        # permition to run bin
        os.system("chmod +700 ./bin_data")
    
    # open file with path to all samples ------------------------------------------
    
    with open("./list_names.txt", "r") as f:
        # create file_content like list elements
        file_content = f.read().split()
        
    
    # path with all combinations of parameters for each sample
    path_base = ["../Results" + second_bar(os.path.dirname(os.path.dirname(file))) for file in file_content]
    
    # all single combinations parameters
    unique_path = list(set([file_base + "/Configs" for file_base in path_base]))
    # auxiliary variable to check if directories Configs are empty
    count = 0
    
    for unique_directories in unique_path:
        all_files = glob.glob(os.path.join(unique_directories,"*.dat"))
        if(len(all_files) == 0):
            print(unique_directories)
            count += 1
        else:
            pass
    
    # if exist Configs folders empty, run script
    if(count != 0):
        # run code for each sample_data
        for file in file_content:
            os.system(f"{exec} {file}")
    else:
        print("all data binning")