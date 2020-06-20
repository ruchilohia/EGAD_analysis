from rank import rank
import numpy as np
import pandas as pd
import bottleneck
import gc
from scipy import sparse
from sklearn.preprocessing import StandardScaler

def compute_sparse_correlation_matrix(A):
    scaler = StandardScaler(with_mean=False)
    scaled_A = scaler.fit_transform(A)  # Assuming A is a CSR or CSC matrix
    corr_matrix = (1/scaled_A.shape[0]) * (scaled_A.T @ scaled_A)
    
    del scaled_A, scaler
    gc.collect()

    if sparse.issparse(corr_matrix):
        return corr_matrix.todense() #Want to return the matrix as dense
    else:
        return corr_matrix

def create_nw(data, replace_nans):

    """Compute Co-expression network from the data
    
    Core network building function. We always run with replace_nans = True
    Slicing single cell data will reguarly produce genes with no counts.
    And any correlation with a vector of all 0s is Nan.
    
    Arguments:
        data {np.array} -- Array of float values, or sparse matrix in shape of genes x cells
        replace_nans {bool} -- Flag for whether to replace Nans in network
    
    Returns:
        np.array -- ranked co-expression matrix of genes x genes 
    """
    nw = compute_sparse_correlation_matrix(data.T)
    np.fill_diagonal(nw, 1)
    nw2 = rank(nw)
    if replace_nans:
        nw2[np.isnan(nw2)] = bottleneck.nanmean(nw2)
    del nw 
    gc.collect()
    return nw2


def nw_aggregation(nw_paths, genes, file_key='nw'):
    """Function for aggregating co-expression networks
    
    Takes a list of paths to HDF5 files and reads in networks,
    avearges them and then re-ranks.

    Each HDF5 needs to be in the Pytable in the fixed format with
    the network stored under the key listed in the keyword argument file_key
    
    Arguments:
        nw_paths {list} -- list of strings or paths to HDF5 files
        genes {np.array} -- numpy array of genes for network
    
    Keyword Arguments:
        file_key {str} --  key in HDF5 network is stored under (default: {'nw'})
    
    Returns:
        pd.DataFrame -- Aggregate Network
    """

    agg_nw = np.zeros([genes.shape[0], genes.shape[0]])
    for nw_path in nw_paths:
        nw = pd.read_hdf(nw_path,file_key)
        fill = bottleneck.nanmean(nw.values,axis=None)
        agg_nw +=nw.loc[genes,genes].fillna(fill).values
        del nw
        gc.collect()

    return pd.DataFrame(rank(agg_nw),index=genes, columns=genes)
