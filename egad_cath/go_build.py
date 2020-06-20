import pandas as pd
new = pd.read_csv('pdb_chain_go.tsv',sep='\t', header=None)

#new=data['PDB,CHAIN,SP_PRIMARY,WITH_STRING,EVIDENCE,GO_ID'].str.split(",", n = 5, expand = True)
new_final=new.drop([2,2,3,4], axis=1) 

protein_list=new[0].to_list()
protein_list_clean=list(set(protein_list)) 
print len(protein_list_clean)
go_list=new[5].to_list()
go_list_clean=list(set(go_list)) 


df_mini = pd.DataFrame(0, index =protein_list_clean, 
                                               columns =go_list_clean)


#from tqdm import tqdm
import numpy as np
dic={}
for i in go_list_clean:
        dic[i]=new[0][new[5]==i]
#data = list(dic.items())
#an_array = np.array(data)

df_mini = pd.DataFrame(0, index =protein_list_clean, 
                                               columns =go_list_clean) 

for i in go_list_clean:
    a=dic[i].to_list()
    for j in a:
        df_mini.at[j,i]=1
#test
#df_mini[df_mini['GO:0032868']==1]
#df_mini
#dic
print df_mini
