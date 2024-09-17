from src.multithread_functions import *


m_relx = [1, 3, 5, 7, 9]*4
j_min_1 = [format(0.98e-03, '.2e') for i in range(5)]
j_min_2 = [format(1.00e-03, '.2e') for i in range(5)]
j_min_3 = [format(1.02e-04, '.2e') for i in range(5)]
j_min_4 = [format(1.04e-04, '.2e') for i in range(5)]
j_min = j_min_1 + j_min_2 + j_min_3 + j_min_4
h_min = j_min

m_teq = [250 for i in range(4*len(j_min_1))]

filename = ["sampleCarmona" for i in range(4*len(j_min))]

method = ['false' for i in range(4*len(j_min))]

multithread_local(filename, j_min, h_min, m_teq, m_relx, method)
#minimum_val()