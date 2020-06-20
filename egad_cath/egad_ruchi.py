import h5py
import gc
import numpy as np
import pandas as pd
from scipy import stats, sparse
import bottleneck
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

from scipy.sparse import csr_matrix

from scipy import sparse

def run_egad(go, nw, **kwargs):
    """EGAD running function
    
    Wrapper to lower level functions for EGAD

    EGAD measures modularity of gene lists in co-expression networks. 

    This was translated from the MATLAB version, which does tiled Cross Validation
    
    The useful kwargs are:
    int - nFold : Number of CV folds to do, default is 3, 
    int - {min,max}_count : limits for number of terms in each gene list, these are exclusive values


    Arguments:
        go {pd.DataFrame} -- dataframe of genes x terms of values [0,1], where 1 is included in gene lists
        nw {pd.DataFrame} -- dataframe of co-expression network, genes x genes
        **kwargs 
    
    Returns:
        pd.DataFrame -- dataframe of terms x metrics where the metrics are 
        ['AUC', 'AVG_NODE_DEGREE', 'DEGREE_NULL_AUC', 'P_Value']
    """
    assert nw.shape[0] == nw.shape[1] , 'Network is not square'
    #print(nw.index)
    #nw.columns = nw.columns.astype(int)
    #print(nw.columns.astype(int))
    assert np.all(nw.index == nw.columns) , 'Network index and columns are not in the same order'

    nw_mask = nw.isna().sum(axis=1) != nw.shape[1]
    nw = nw.loc[nw_mask, nw_mask].astype(float)
    np.fill_diagonal(nw.values, 1)
    return _runNV(go, nw, **kwargs)


def _runNV(go, nw, nFold=3, min_count=0, max_count=100000):

    #Make sure genes are same in go and nw
    genes_intersect = go.index.intersection(nw.index)

    go = go.loc[genes_intersect, :]
    nw = nw.loc[genes_intersect, genes_intersect]

    #Make sure there aren't duplicates
    duplicates = nw.index.duplicated(keep='first')
    nw = nw.loc[~duplicates, ~duplicates]

    go = go.loc[:, (go.sum(axis=0) > min_count) & (go.sum(axis=0) < max_count)]
    go = go.loc[~go.index.duplicated(keep='first'), :]
    #print(go)

    roc = _new_egad(go.values, nw.values, nFold)

    col_names = ['AUC', 'AVG_NODE_DEGREE', 'DEGREE_NULL_AUC', 'P_Value']
    #Put output in dataframe
    return pd.DataFrame(dict(zip(col_names, roc)), index=go.columns)


def _new_egad(go, nw, nFold):

    #Build Cross validated Positive
    x, y = np.where(go)
    #print(x, y)
    cvgo = {}
    for i in np.arange(nFold):
        a = x[i::nFold]
        #print(a)
        b = y[i::nFold]
        dat = np.ones_like(a)
        mask = sparse.coo_matrix((dat, (a, b)), shape=go.shape)
        cvgo[i] = go - mask.toarray()
        
    CVgo = np.concatenate(list(cvgo.values()), axis=1)
    #print(CVgo)

    sumin = np.matmul(nw.T, CVgo)

    degree = np.sum(nw, axis=0)
    #print(degree)
    #print(degree[:, None])

    predicts = sumin / degree[:, None]
    #print(predicts)

    np.place(predicts, CVgo > 0, np.nan)

    #print(predicts)

    #Calculate ranks of positives
    rank_abs = lambda x: stats.rankdata(np.abs(x))
    predicts2 = np.apply_along_axis(rank_abs, 0, predicts)
    #print(predicts2)

    #Masking Nans that were ranked (how tiedrank works in matlab)
    predicts2[np.isnan(predicts)] = np.nan
    #print(predicts2)

    filtering = np.tile(go, nFold)
    #print(filtering)

    #negatives :filtering == 0
    #Sets Ranks of negatives to 0
    np.place(predicts2, filtering == 0, 0)

    #Sum of ranks for each prediction
    p = bottleneck.nansum(predicts2, axis=0)

    #Number of predictions
    #Number of 1's masked for each GO term for each CV
    n_p = np.sum(filtering, axis=0) - np.sum(CVgo, axis=0)

    #Number of negatives
    #Number of GO terms - number of postiive
    n_n = filtering.shape[0] - np.sum(filtering, axis=0)

    roc = (p / n_p - (n_p + 1) / 2) / n_n
    U = roc * n_p * n_n
    Z = (np.abs(U - (n_p * n_n / 2))) / np.sqrt(n_p * n_n *
                                                (n_p + n_n + 1) / 12)
    roc = roc.reshape(nFold, go.shape[1])
    Z = Z.reshape(nFold, go.shape[1])
    #Stouffer Z method
    Z = bottleneck.nansum(Z, axis=0) / np.sqrt(nFold)
    #Calc ROC of Neighbor Voting
    roc = bottleneck.nanmean(roc, axis=0)
    P = stats.norm.sf(Z)

    #Average degree for nodes in each go term
    avg_degree = degree.dot(go) / np.sum(go, axis=0)

    #Calc null auc for degree
    ranks = np.tile(stats.rankdata(degree), (go.shape[1], 1)).T

    np.place(ranks, go == 0, 0)

    n_p = bottleneck.nansum(go, axis=0)
    nn = go.shape[0] - n_p
    p = bottleneck.nansum(ranks, axis=0)

    roc_null = (p / n_p - ((n_p + 1) / 2)) / nn
    print roc 
    return roc, avg_degree, roc_null, P

#go = pd.read_csv("/home/lohia/egad_trial/annotations")
#go = pd.read_hdf('/data/bharris/GO_data/go_mouse_nw.hdf5', 'go')

with h5py.File("/data/johlee/CoCoCoNet/gene2go/yeast_gene2go.hdf5", "r") as f:
    print f.keys()   # works like a dict

    a_group_key = list(f.keys())[2]
    print list(f[a_group_key])
#print pd.read_hdf('/data/johlee/CoCoCoNet/gene2go/yeast_gene2go.hdf5', ['GO', 'gene']) 
#print go
#nw = pd.read_csv("/home/lohia/egad_trial/gene_network")
#nw = pd.read_hdf('/data/bharris/biccn_coexpression/data/bulk_rna/aggregate_bulk_nw.hdf5', 'nw')
#print nw

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
data['cath']=new[0]
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
print(data.shape)
#data.drop_duplicates(subset='etc.', keep='first', inplace=True)
#print data

protein_list=data['pdb_chain_domain'].to_list()
num_protein = len(protein_list)
print len(protein_list)
pre_values=data.iloc[:,1].values
#print pre_values.shape
#p_values = np.reshape(pre_values,(num_protein,1))
#data_matrix =  np.tile(pre_values,(num_protein,1))
data_matrix = np.zeros((num_protein,num_protein), dtype='int8')
#data_matrix = data_matrix.astype('int8') 
#data_matrix = (sparse.lil_matrix(data_matrix))
#print data_matrix
#a1,b1,c1,d1,e1,f1,g1=  np.array_split(data_matrix,7)
#p1,p2,p3,p4,p5,p6,p7 = np.array_split(pre_values,7)
#print p1
#output_matrix = data_matrix / pre_values[:, None]
row_idx =0
for cluster in pre_values:
       data_matrix[row_idx] = pre_values / cluster
       data_matrix[row_idx] = np.where(data_matrix[row_idx]==1, 1,0)
       row_idx = row_idx +1
print(data_matrix)

output_matrix = data_matrix
#output_matrix = data_matrix.multiply(rD_sp)
#print(output_matrix)
print output_matrix
#o_p_2 = b1 / p2[:,None] 
#print output_matrix.tolist()
#output_matrix = np.concatenate((output_matrix_1, o_p_2))
#nw = np.where(output_matrix==1, 1,0)
#print nw
nw = pd.DataFrame(data=output_matrix,index=protein_list,columns=protein_list)
print nw

#nw.to_csv('nw_C.csv')
nw.to_hdf('nw_C.h5', key='nw', mode='w')
gc.collect()
#new=data['PDB,CHAIN,SP_PRIMARY,WITH_STRING,EVIDENCE,GO_ID'].str.split(",", n = 5, expand = True)
#new=data['PDB,CHAIN,SP_PRIMARY,WITH_STRING,EVIDENCE,GO_ID'].str.split(",", n = 5, expand = True)
new = pd.read_csv('pdb_chain_go.tsv',sep='\t', header=None, nrows=100)


new_final=new.drop([2,2,3,4], axis=1) 

protein_list=new[0].to_list()
protein_list_clean=list(set(protein_list)) 
print len(protein_list_clean)
go_list=new[5].to_list()
go_list_clean=list(set(go_list)) 


df_mini = pd.DataFrame(0, index =protein_list_clean, 
                                               columns =go_list_clean)



dic={}
for i in go_list_clean:
        dic[i]=new[0][new[5]==i]
df_mini = pd.DataFrame(0, index =protein_list_clean, 
                                               columns =go_list_clean) 

for i in go_list_clean:
    a=dic[i].to_list()
    for j in a:
        df_mini.at[j,i]=1
#print df_mini
#df_mini.to_hdf('go.h5', key='df', mode='w')
df_mini = pd.read_hdf('go.h5', 'df')
print df_mini.head()
df = run_egad(df_mini, nw)
#print df[['AUC']]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 4), sharey=True)
ax = df.plot.scatter(x='AUC',y='DEGREE_NULL_AUC')
print df.mean()
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.plot([0, 1], [0, 1], transform=ax.transAxes, c='black')
plt.axvline(x=df['AUC'].mean(),c='black',ls='--')
plt.axhline(y=df['DEGREE_NULL_AUC'].mean(), c='black', ls='--')
plt.savefig('C.pdf', bbox_inches='tight', dpi=100)



#df.to_csv('C.csv')
