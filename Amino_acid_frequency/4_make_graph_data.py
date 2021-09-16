#%%
import pandas as pd
import glob

protein_names = ['N', 'NS9c', 'NSP12', 'NSP2'] #include clades that have differences in the most frequent amino acids
clade_names = ['G', 'GH', 'GR', 'GV', 'L', 'O', 'S', 'V'] 

gisaid_file = glob.glob("path to the directory that have the files of amino acid frequencies")

for clade_folder, clade_name in zip(gisaid_file, clade_names):
    protein_file = glob.glob('{}/*'.format(clade_folder))
    for graph_file in protein_file:
        if protein_file[2] == graph_file or protein_file[9] == graph_file or protein_file[13] == graph_file or protein_file[18] == graph_file:
            if protein_file[2] == graph_file: #if the file is the file of N protein
                protein_name = protein_names[0]
            elif protein_file[9] == graph_file: #if the file is the file of NS9c protein
                protein_name = protein_names[1]
            elif protein_file[13] == graph_file:
                protein_name = protein_names[2]
            elif protein_file[18] == graph_file: #if the file is the file of NSP2 protein
                protein_name = protein_names[3]  

            df = pd.read_csv(graph_file)

            df2 = df.iloc[:, 1::2]
            df3 = df.iloc[:, 2::2]

            dic = {}
            lis = []
            len_columns = (len(df2.columns))
            for (index2, items2), (index3, items3) in zip(df2.iterrows(), df3.iterrows()):
                for i in range(len_columns):
                    dic[items2[i]] = items3[i]
                lis.append(dic)
                dic = {}

            df_fin = pd.DataFrame(lis)
            df_fin.insert(0, 'position', df_fin.index+1)
            df_fin.dropna(how ='all', axis=1, inplace=True) 
            df_fin.to_csv("name and path of the csv file", index=False)


# %%
