from src.multithread_functions import *

m_relx =  [i for i in range(3, 10, 2) for _ in range(5)]
m_teq = [250, 300, 350, 400, 450] * 4
seed = [i for i in range(20,40)]
filename = ["data_synteticN40" for i in range(len(m_relx))]
method = ["metropolis" for i in range(len(m_relx))]
j_min = [format(7.89e-04, '.2e') for i in range(len(m_relx))]
h_min = [format(4.85e-04, '.2e') for i in range(len(m_relx))]
comb_parms = [(i, j) for i, j in zip(m_teq, m_relx)]


multithread_local(filename, j_min, h_min, m_teq, m_relx, seed, method)



