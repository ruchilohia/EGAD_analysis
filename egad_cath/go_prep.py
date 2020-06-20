import pandas as pd
import pickle
import numpy as np
#import h5py

data = pd.read_csv('goa_pdb.gaf',sep='\t', header=1,skiprows=range(1, 16))

data_subset=data.loc[:, ['101M_A','GO:0005344','taxon:9755']]
data_subset.columns = ['PDB_ID','GO_ID','TAXON']
#taxonomy for human 9606
data_subset_human = data_subset[data_subset['TAXON']=='taxon:9606']  
data_subset_human = data_subset_human.drop_duplicates(subset=['PDB_ID','GO_ID'])
data_subset_human

pickle_in = open('gotermindex.pickle','rb')
pickle_in2 = open('altiddict.pickle','rb')
gotermidx_dict = pickle.load(pickle_in)
altid_dict = pickle.load(pickle_in2)
ontology_np = np.load('ontology.npy')


list_of_pdb=data_subset_human['PDB_ID'].tolist()
list_of_pdb_dedupe = list( dict.fromkeys(list_of_pdb) )
i = 0
pdbidx = {}
for pdb in list_of_pdb_dedupe:
    pdbidx[pdb] = i
    i = i + 1


shape = (len(list_of_pdb_dedupe), len(ontology_np[1]))
pdbmatrix = np.zeros( shape, dtype=bool )
pdbmatrix.shape


pdb_go_dict={}
ct=0
for i in list_of_pdb_dedupe:
    minidf=data_subset_human[data_subset_human['PDB_ID']==i]
    list_of_gos=minidf['GO_ID'].tolist()
    #print(len(list_of_gos))
    list_of_gos_idx=[]
    for j in list_of_gos:
        try:
            list_of_gos_idx.append(gotermidx_dict[j])  # use atl id if go id not found in go term id
        except:
            ALT_PARENT_GO_ID=altid_dict[j]
            list_of_gos_idx.append(gotermidx_dict[ALT_PARENT_GO_ID]) 
            #list_of_gos_idx.append() 
    pdb_go_dict[ct]=list_of_gos_idx
    ct=ct+1


for i in pdb_go_dict.keys():
     for j in pdb_go_dict[i]:
            pdbmatrix[i][j] = True
print(pdbmatrix)
go_id = gotermidx_dict.keys()
print(go_id)
np.save('go_id.npy' , list(go_id))
np.save('go_pdb.npy', pdbmatrix)
np.save('pdb_index.npy', list_of_pdb_dedupe)
#np.save('go_id', go_id)
#f = open('pdb_index', 'wb')   
#pickle.dump(list_of_pdb_dedupe, f, protocol=2 )
#f.close()
#df_mini = pd.DataFrame(pdbmatrix, index =list_of_pdb_dedupe,  columns = go_id) 
#df_mini.to_hdf('go_prop.h5', key='df', mode='w')
#print(df_mini)
