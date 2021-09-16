import pandas as pd
import os

#exclude unnecessary columns from meta data of GISAID clade
data_frame = pd.read_csv("path of metadata downloaded from GISAID", usecols=[2,13,19], index_col=False)
data_frame.set_index("gisaid_epi_isl")
data_frame.to_csv("csv file of metadata", index=False)

#exclude unnecessary columns from fasta file data of GISAID clade
df = pd.read_csv("path of FASTA file downloaded from GISAID", header=None, lineterminator = '>', sep='|', usecols = [0,3,10])
df[10] = df[10].str.split('\n', expand = True)[1]
df.columns = ['protein_name', 'gisaid_epi_isl', 'sequence']
df.to_csv("csv file of sequence data", index=False)

#create joint file of metadata and sequence data
df = pd.read_csv("csv file of sequence data")
df2 = pd.read_csv("csv file of metadata")
data_frame = pd.merge(df, df2, on="gisaid_epi_isl", how='right')
data_frame.to_csv("csv file of sequence data including sequence ID and clade data", index=False)

#create FASTA file for each clade and each protein
clade_name = ["G", "GR", "GH", "GV", "S", "O", "L", "V"]
for clade in clade_name:
    os.mkdir("directory for each clade")
    for protein in data_frame['protein_name'][0:27]:
        #just include sequences longer than 29000 bases
        df2 = data_frame.loc[(data_frame['GISAID_clade'] == clade) & (data_frame['protein_name'] == protein) & (data_frame['length'] > 29000)]
        new_df2 = df2.drop_duplicates()
        new_df3 = new_df2.reset_index(drop=True)
        new_df3 = new_df3.drop(columns=['length'])
        new_df3 = new_df3.iloc[:, [0,3,1,2]]

        new_df3['protein_name'] = '>' + new_df3['protein_name']
        new_df3['sequence'] = '\n' + new_df3['sequence']

        new_df3.to_csv("path of FASTA file", index=False, header=False, sep='|')

        #remove doublequotes
        myfile = open("path of FASTA file", 'r')
        data = myfile.read()
        data = str(data).replace('"', "")
        data = str(data).replace('*', "")
        with open("path of edited FASTA file", 'w') as f:
            f.write(data)
        myfile.close()