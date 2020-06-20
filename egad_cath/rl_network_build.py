import pandas as pd
import numpy as np
data = pd.read_csv('cath-b-newest-all.txt', sep=" ", header=None, nrows=30)

data.columns = ["a", "b", "c", "etc."]
#print data
data['a']=data['a'].str[:4] # removing extra characters in protein name
new=data['c'].str.split(".", n = 3, expand = True)
data['c']=new[0]
data['c'] = data['c'].astype(int)
#data[5]=new[1]
#data[6]=new[2]
#data[7]=new[3]

#protein_list_clean = list(set(protein_list)) 
#protein_list_mini=protein_list_clean[:100]

del data['b']
del data['etc.']
#del data[3]
#print data
data.drop_duplicates(subset='a', keep='first', inplace=True)
protein_list=data['a'].to_list()
num_protein = len(protein_list)
print len(protein_list)
pre_values=data.iloc[:,1].values
#print pre_values.shape
#p_values = np.reshape(pre_values,(num_protein,1))
data_matrix =  np.tile(pre_values,(num_protein,1))
#print data_matrix
#a1,b1,c1,d1,e1,f1,g1=  np.array_split(data_matrix,7)
#p1,p2,p3,p4,p5,p6,p7 = np.array_split(pre_values,7)
#print p1
output_matrix = data_matrix / pre_values[:, None]
#o_p_2 = b1 / p2[:,None] 
#print output_matrix.tolist()
#output_matrix = np.concatenate((output_matrix_1, o_p_2))
nw = np.where(output_matrix==1, 1,0)
print nw
df = pd.DataFrame(data=nw,index=protein_list,columns=protein_list)
print df
