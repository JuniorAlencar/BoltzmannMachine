from src.directories_fluids import *
import glob
#copy_directory_structure(src_dir, dest_dir)
#find_dat_files_in_configs(root_dir, output_file)

# Exemplo de uso

# Função para pegar o caminho após a segunda barra
def pegar_apos_segunda_barra(caminho):
    partes = caminho.split('/')
    return '/' + '/'.join(partes[3:]) if len(partes) > 3 else caminho


if __name__ == "__main__":
    # reproduce structure of folders to flow_spin data
    src_directory = '../Data/BruteData'  # folder with brute data
    dest_directory = '../Results'        # folders to binning data

    copy_directory_structure(src_directory, dest_directory)
    # create file with path to all files .dat to binning ---------------------
    root_dir = "../Data/BruteData"
    output_file = "./list_names.txt"

    find_dat_files_in_configs(root_dir, output_file)
    
    # create executable(if dont exist) to binning data in c++ -------------------------
    exec = "./bin_data"
    if os.path.isfile(exec):
        pass
    else:
        os.system("g++ -std=c++17 -o bin_data ../src/clean_data.cpp")
        # permition to run bin
        os.system("chmod +700 ./bin_data")
    
    # open file with path to all samples ------------------------------------------
    
    with open("list_names.txt", "r") as f:
        # create file_content like list elements
        file_content = f.read().split()
    # path with all combinations of parameters for each sample
    path_base = [pegar_apos_segunda_barra(os.path.dirname(os.path.dirname(file))) for file in file_content]
    print(path_base)
    # all single combinations parameters
    unique_path = list(set([file_base + "/Configs" for file_base in path_base]))
    # auxiliary variable to check if directories Configs are empty
    count = 0
    
    for unique_directories in unique_path:
        all_files = glob.glob(os.path.join(unique_directories))
        print(all_files)
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
        
        

        

    
