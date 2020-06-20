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