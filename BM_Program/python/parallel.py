from src.multithread_functions import *

j_min_1 = [format(1.03e-03, '.2e'), format(1.00e-03, '.2e'), format(9.98e-04, '.2e'), format(9.96e-04, '.2e'),format(9.94e-04, '.2e')]
h_min_1 = [format(1.03e-03, '.2e'), format(1.00e-03, '.2e'), format(9.98e-04, '.2e'), format(9.96e-04, '.2e'),format(9.94e-04, '.2e')]

j_min = j_min_1*4
h_min = h_min_1*4

m_teq_1 = [350 for i in range(len(j_min_1))]
m_teq_2 = [300 for i in range(len(j_min_1))]
m_teq_3 = [250 for i in range(len(j_min_1))]
m_teq_4 = [200 for i in range(len(j_min_1))]

m_teq = m_teq_1 + m_teq_2 + m_teq_3 + m_teq_4

filename = ["sampleCarmona" for i in range(4*len(j_min))]
m_relx = [2 for i in range(4*len(j_min))]
method = ['false' for i in range(4*len(j_min))]

multithread_local(filename, j_min, h_min, m_teq, m_relx, method)