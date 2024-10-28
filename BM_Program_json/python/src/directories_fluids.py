import os
import shutil

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
                        # Armazenar o caminho completo do arquivo
                        dat_files.append(os.path.join(root, file))

        # Escrever o caminho dos arquivos .dat em um arquivo de texto
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
            # Criar o caminho correspondente no diretório de destino
            relative_path = os.path.relpath(root, src_dir)  # Caminho relativo da estrutura
            dest_path = os.path.join(dest_dir, relative_path)  # Caminho de destino

            # Criar o diretório no destino
            os.makedirs(dest_path, exist_ok=True)
            print(f"Diretório criado: {dest_path}")
