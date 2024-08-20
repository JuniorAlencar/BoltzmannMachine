import os
import glob

# sample data: filename sample without extesion sampleN20.dat -> sampleN20
# h_min: minimum value h to MC converge
# j_min: minimum value J to MC converge
# exact_solutions: if True use, if False use Metropolis
# eta_J and eta_h: process update rate to J and h respectively
# iter_max: number maximum number of interactions in MC
# n_rep: ?
# multiply_teq: t_eq = nspins * multiply_teq
def gen_json_input(sample_data, h_min, j_min, exact_solutions=True ,eta_J=0.05, eta_h=0.03, iter_max=150000, n_rep=40, multiply_teq=150):
    
    # Create folder inputs
    newpath = f"../inputs/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    input_data = "../Data/" + sample_data + ".dat"

    with open(input_data, 'r') as file:
        # Read the first line
        first_line = file.readline().strip()

        # Split the line into columns based on whitespace
        columns = first_line.split(',')

        # Count the number of columns
        nspins = len(columns)

    t_eq = multiply_teq * nspins
    t_rep = 2 * nspins
    t_meas = int(nspins * 6000 * t_rep / n_rep)

    if(exact_solutions == False):
        method = "metropolis"
        exact_solutions = 'false'
        folder_results = "../Results_Metropolis/"
    elif(exact_solutions == True):
        method = "exact"
        exact_solutions = 'true'
        folder_results = "../Results/"

    filename = f"input_j_min_{j_min}_h_min_{h_min}_t_eq_{t_eq}_{method}.json"

    a= "{ \n"
    b = f'"run_name":"{sample_data}",\n'
    c = f'"input_data":"../Data/{sample_data}.dat",\n'
    d = f'"input_init_guess":' + f'{folder_results}/{sample_data},\n'
    d = f'"input_init_guess":"../Results/{sample_data}", \n'
    e = f'"output_err": "../Results/err_{sample_data}_j_min_{j_min}_h_min_{h_min}_t_eq_{t_eq}_{method}.txt", \n'
    f = f'"output_data": "../Results/{sample_data}_j_min_{j_min}_h_min_{h_min}_t_eq_{t_eq}_{method}.json", \n'
    g = f'"use_exact": {exact_solutions},\n\n'
    j = '"relax": {\n'
    h = '\t\t\t"iter_display": 100, \n'
    l = f'\t\t\t"iter_max": {iter_max}, \n'
    m = f'\t\t\t"eta_J": {eta_J}, \n'
    n = f'\t\t\t"eta_h": {eta_h}, \n'
    o = f'\t\t\t"min_err_S":{h_min},\n'
    p = f'\t\t\t"min_err_SS":{j_min}\n'
    q = "\t\t\t\t},\n\n"
    s = '"mc": {\n'
    u = f'\t\t\t"n_rep": {n_rep},\n'
    v = f'\t\t\t"t_eq": {t_eq},\n'
    w = f'\t\t\t"t_meas": {t_meas}\n'
    x = "\t\t\t}\n"
    y = "}"
    lst = [a+b+c+d+e+f+g+j+h+l+m+n+o+p+q+s+u+v+w+x+y]
    x = open(newpath + filename, "w")
    for k in lst:
        x.write(k)
    x.close()

def gen_shell_multithread():
    filename = f"multithread_pc.sh"

    a = "#!/bin/bash\n\n"

    b = "# Define uma função que contêm o código para rodar em paralelo\n"

    c = "run_code() {\n\t"
    d = f"time ../bin/bmc ../inputs/$1\n"
    e = "}\n"
    f = "# Exportar a função usando o módulo Parallel\n"
    g = "export -f run_code\n\n"

    path_d = f"../inputs"
    all_files = glob.glob(os.path.join(path_d,"*.json"))
    list_of_arguments = [V[2] for V in os.walk(path_d)][0]
    list_of_arguments = str(list_of_arguments)
    list_of_arguments = list_of_arguments.replace(',', '')

    h = f"arguments=(" 
    i = list_of_arguments[1:-1] + ")\n"
    j = "parallel run_code :::\t" +  """ "${arguments[@]}"  """ "\n\t"
    list_for_loop = [a,b,c,d,e,g,h,i,j]
    l = open("../Scripts/" + filename, "w") # argument w: write if don't exist file
    for k in list_for_loop:
        l.write(k)
    l.close()

def permission_run():
    filename = f"../multithread_pc.sh"
    os.system(f"chmod 700 ../Scripts/" + filename)