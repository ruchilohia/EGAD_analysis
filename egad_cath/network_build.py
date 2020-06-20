import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csr_matrix

def build_netowk(id_val=0):

    data = pd.read_csv('cath-b-newest-all.txt', sep=" ", header=None) 

    data.columns = ["pdb_chain_domain", "version", "cath", "residues"]
    del data['version']
    del data['residues']
    #print data
    data['pdb']=data['pdb_chain_domain'].str[:4] # removing extra characters in protein name
    data['pdb'] = data['pdb'].apply(lambda x: x.upper())
    data['pdb_chain'] = data['pdb'] + '_' + data['pdb_chain_domain'].str[4:5]
    #del data['pdb_chain_domain']   

    #del data['pdb']
    new=data['cath'].str.split(".", n = 3, expand = True)  # expand the cath domains
    #print(data.to_string())
    data['cath']=new[id_val]
    #data['etc.']=new[3] #just takes same homology family into account
    data['cath'] = data['cath'].astype(float)
    #data[5]=new[1]
    #data[6]=new[2]
    #data[7]=new[3] 

    #protein_list_clean = list(set(protein_list)) 
    #protein_list_mini=protein_list_clean[:100] 
    

    #del data[3]
    #print(data.shape)
    data.drop_duplicates(subset=['pdb', 'cath'], keep='first', inplace=True)
    del data['pdb']
    #print(data.shape)
    #data.drop_duplicates(subset='etc.', keep='first', inplace=True)
    #print data 

    protein_list=data['pdb_chain'].to_list()
    num_protein = len(protein_list)
    #print len(protein_list)
    pre_values=data.iloc[:,1].values
    pre_values = pre_values.astype('int8')
    #print pre_values.shape
    #p_values = np.reshape(pre_values,(num_protein,1))
    #data_matrix =  np.tile(pre_values,(num_protein,1))
    data_matrix = np.zeros((num_protein,num_protein), dtype='int8')
    #print(data_matrix.dtype)   
    

    row_idx =0
    for cluster in pre_values:
           #col_idx = row_idx +1
           #data_matrix[row_idx][0 :col_idx]   = 5
           #data_matrix[row_idx][0 :col_idx] = pre_values[0 :col_idx] / cluster
           data_matrix[row_idx] = np.where(pre_values == cluster, 1,0)
           #data_matrix= np.maximum( data_matrix, data_matrix.transpose() )
           row_idx = row_idx +1
    #data_matrix = data_matrix.astype('float64')
    output_matrix = data_matrix
    #print output_matrix
    #nw = pd.DataFrame(data=output_matrix,index=protein_list,columns=protein_list)
    np.save("nw.npy", output_matrix)    
#nw.to_hdf('nw_%s.h5' %id_val, key='nw', mode='w')
build_netowk(0)
