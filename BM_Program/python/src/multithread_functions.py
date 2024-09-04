import os

#filename: filename of data (without extension)
#j_min: minimum value Jij to MC converge (in scientific notation)
#h_min: minimum value hi to MC converge (in scientific notation)
#m_relx: multiply relx (relx = nspins*m_relx)
def multithread_local(filename, j_min, h_min, m_teq, m_relx, method):
    # file name with inputs
    file_inputs = "input_multithread.txt"
    folder = "../Scripts/"
    
    # write .txt with inputs to parallel
    # Ensure single space between values
    lst = [' '.join(map(str, (i, j, k, l, m, n))) for i, j, k, l, m, n in zip(filename, j_min, h_min, m_teq, m_relx, method)]
    
    # writing inputs
    with open(folder + file_inputs, "w") as file:
        for line in lst:
            file.write(line + "\n")
    
    # generate the outputs to run code
    
    # parallel script
    with open(folder + "parallel.sh", 'w') as a:
        # generate the outputs to run code
        a.write("#!/bin/bash\n\n")
        a.write("# Generate the executables\n")
        a.write("./compile.sh\n\n")
        a.write("# Create folders to Results if don't exist\n")
        a.write('if [ ! -d "../Results" ] || [ ! -d "../Results_Metropolis" ]; then\n')
        a.write(f"\t./CreateFolders\n")
        a.write("fi\n\n")
        a.write(f"# Run in parallel with arguments in {file_inputs}\n")
        a.write(f"parallel --colsep ' ' ./main_single.sh {{1}} {{2}} {{3}} {{4}} {{5}} {{6}} < {file_inputs}")
    
    # permission to run parallel script
    os.system("chmod +x ../Scripts/parallel.sh")
    a.close
