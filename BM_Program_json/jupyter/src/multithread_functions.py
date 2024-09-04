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
    os.system("../Scripts/compile.sh")

    # parallel script
    with open(folder + "parallel.sh", 'w') as a:
        a.write("./compile.sh\n\n")
        a.write(f"parallel ./main_single.sh {{1}} {{2}} {{3}} {{4}} {{5}} {{6}} < {file_inputs}")
    
    # permission to run parallel script
    os.system("chmod +x ../Scripts/parallel.sh")

