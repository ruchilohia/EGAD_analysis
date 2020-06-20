data = pd.read_csv("cath-b-newest-all.txt", sep=" ", header=None)

data.columns = ["pdb_chain_domain", "version", "cath", "residues"]
del data["version"]
del data["residues"]
# print data
data["pdb"] = data["pdb_chain_domain"].str[
    :4
]  # removing extra characters in protein name
data["pdb"] = data["pdb"].apply(lambda x: x.upper())
data["pdb_chain"] = data["pdb"] + "_" + data["pdb_chain_domain"].str[4:5]
# del data['pdb_chain_domain']

# del data['pdb']
new = data["cath"].str.split(".", n=3, expand=True)  # expand the cath domains
# print(data.to_string())
data["cath"] = new[0]
# data['etc.']=new[3] #just takes same homology family into account
data["cath"] = data["cath"].astype(float)
# data[5]=new[1]
# data[6]=new[2]
# data[7]=new[3]

# protein_list_clean = list(set(protein_list))
# protein_list_mini=protein_list_clean[:100]


# del data[3]
# print(data.shape)
data.drop_duplicates(subset=["pdb", "cath"], keep="first", inplace=True)
del data["pdb"]
print (data.shape)
# data.drop_duplicates(subset='etc.', keep='first', inplace=True)
# print data

protein_list = data["pdb_chain_domain"].to_list()
num_protein = len(protein_list)
print len(protein_list)
pre_values = data.iloc[:, 1].values
# print pre_values.shape
# p_values = np.reshape(pre_values,(num_protein,1))
# data_matrix =  np.tile(pre_values,(num_protein,1))
data_matrix = np.zeros((num_protein, num_protein), dtype="int8")
# data_matrix = data_matrix.astype('int8')
# data_matrix = (sparse.lil_matrix(data_matrix))
# print data_matrix
# a1,b1,c1,d1,e1,f1,g1=  np.array_split(data_matrix,7)
# p1,p2,p3,p4,p5,p6,p7 = np.array_split(pre_values,7)
# print p1
# output_matrix = data_matrix / pre_values[:, None]
row_idx = 0
for cluster in pre_values:
    data_matrix[row_idx] = pre_values / cluster
    data_matrix[row_idx] = np.where(data_matrix[row_idx] == 1, 1, 0)
    row_idx = row_idx + 1
print (data_matrix)

output_matrix = data_matrix
# output_matrix = data_matrix.multiply(rD_sp)
# print(output_matrix)
print output_matrix
# o_p_2 = b1 / p2[:,None]
# print output_matrix.tolist()
# output_matrix = np.concatenate((output_matrix_1, o_p_2))
# nw = np.where(output_matrix==1, 1,0)
# print nw
nw = pd.DataFrame(data=output_matrix, index=protein_list, columns=protein_list)
print nw

# nw.to_csv('nw_C.csv')
nw.to_hdf("nw_C.h5", key="nw", mode="w")
