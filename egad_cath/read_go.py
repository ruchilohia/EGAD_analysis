import pickle
import numpy as np


pickle_in = open('gotermindex.pickle','rb')

gotermidx_dict = pickle.load(pickle_in)
go_id = gotermidx_dict.keys()

go_matrix = np.load('go_pdb.npy')

pickle_in = open('pdb_index','rb')

pdb_id = pickle.load(pickle_in)
print(pdb_id)

