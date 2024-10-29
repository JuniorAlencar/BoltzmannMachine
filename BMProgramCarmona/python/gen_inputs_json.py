from src.gen_inputs_functions import *

# sample data: filename sample without extesion sampleN20.dat -> sampleN20
# h_min: minimum value h to MC converge
# j_min: minimum value J to MC converge
# exact_solutions: if True use, if False use Metropolis
# eta_J and eta_h: process update rate to J and h respectively
# iter_max: number maximum number of interactions in MC
# n_rep: ?
# multiply_teq: t_eq = nspins * multiply_teq

sample_data = "sampleN20"
h_min = [1.00e-04,5.81e-05]
j_min = [1.00e-04,9.85e-06]
solu = [False, True]

clear_input_folder()

for i in range(len(h_min)):
    gen_json_input(sample_data, h_min[i], j_min[i], exact_solutions=solu[i] ,eta_J=0.05, eta_h=0.03, iter_max=150000, n_rep=40, multiply_teq=150)

gen_shell_multithread()
permission_run()