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
h_min = 4.2e-05
j_min = 9.1e-05
h_min = [3.8e-05, 3.6e-05, 3.4e-05, 3.2e-05]
j_min = [8.5e-05, 8.3e-05, 8.1e-05, 7.9e-05]
for i in range(3):
    for j in range(3):
        gen_json_input(sample_data, h_min[i], j_min[j], exact_solutions=False ,eta_J=0.05, eta_h=0.03, iter_max=150000, n_rep=40, multiply_teq=150)

gen_shell_multithread()
permission_run()