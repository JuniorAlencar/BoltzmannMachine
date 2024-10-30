from src.directories_fluids import *

# copy_directory_structure(src_dir, dest_dir)
# find_dat_files_in_configs(root_dir, output_file)
# second_bar(path)
# binning_data()

if __name__ == "__main__":
    # reproduce structure of folders to flow_spin data
    src_directory = '../Data/BruteData'  # folder with brute data
    dest_directory = '../Results'        # folders to binning data
    
    copy_directory_structure(src_directory, dest_directory)
    
    # create file with path to all files .dat to binning ---------------------
    root_dir = "../Data/BruteData"
    output_file = "./list_names.txt"
    
    # find all paths in configs samples and save in list_names.txt file
    find_dat_files_in_configs(root_dir, output_file)
    
    # binning all data
    binning_data()

    # Calculate all experimental means from samples
    calculate_all_means()
    
    
        
        

        

    
