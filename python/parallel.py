from src.multithread_functions import *

m_relx = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
m_teq = [100, 150, 200, 250] * 4
#seed = [1,2,3,4,5,6,7,8,9,10,11,12]
comb_parms = [(i, j) for i, j in zip(m_teq, m_relx)]

filename = "data_synteticN40"
method = "metropolis"


