from src.multithread_functions import *


m_relx = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
j_min = [format(1.60e-03, '.2e')]*12
h_min = [format(1.00e-03, '.2e')]*12
m_teq = [100,150, 200, 250]*4
seed = [1,2,3,4,5,6,7,8,9,10,11,12]

filename = ["data_synteticN30" for i in range(len(m_relx))]

method = ['metropolis' for i in range(len(m_relx))]

#multithread_local(filename, j_min, h_min, m_teq, m_relx, seed, method)
minimum_val("metropolis", seed)