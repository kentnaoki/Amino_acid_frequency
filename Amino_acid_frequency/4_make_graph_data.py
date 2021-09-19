#%%
import pandas as pd
import glob
import re
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='create a graph representing amino acid frequency in GISAID clade')
parser.add_argument('protein_name', help='input protein name')
parser.add_argument('aa_change', help='input amino acid change (e.g. R203K)')
args = parser.parse_args()

#function to create data for graph
def make_graph_data(protein):
    gisaid_file_freq_posi = glob.glob("/mnt/c/Users/naoki/Onedrive/デスクトップ/MSc_research/clade_frequency/Frequency_results/clade_aa_frequency_correct/GISAID/**/*.*", recursive=False) #files must be csv files
    for graph_file in gisaid_file_freq_posi:
        graph_file_protein = re.split('[/_.]', graph_file)[-2] #The file name must be 'CladeName_ProteinName.csv'
        graph_file_clade = re.split('[/_.]', graph_file)[-3]
        if graph_file_protein == protein:
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
            df_fin.to_csv("/mnt/c/Users/naoki/Onedrive/デスクトップ/github_source_code/{}_{}.csv".format(graph_file_clade, graph_file_protein), index=False)

#function to create a graph
def create_graph(protein, aa_change):
    clade_names = []
    position = int(aa_change[1:-1])
    freq_list = []
    gisaid_file_freq_aa = glob.glob("/mnt/c/Users/naoki/Onedrive/デスクトップ//github_source_code/*.csv") #files must be csv files
    for graph_file in gisaid_file_freq_aa:
        graph_file_protein = re.split('[/_.]', graph_file)[-2] #The file name must be 'CladeName_ProteinName.csv'
        graph_file_clade = re.split('[/_.]', graph_file)[-3]
        if graph_file_protein == protein:
            clade_names.append(graph_file_clade)
            df = pd.read_csv(graph_file)
            df2 = df.drop(df.columns[0], axis=1)
            df2 = df2.fillna(0).astype('int64', errors='ignore')
            df2 = df2.astype(float)
        #calculate the ratio of amino acids in the position
            sum_count = 0
            for m in range(len(df2.iloc[0])):
                sum_count += df2.iloc[0][m]
            for j in range(len(df2.index)):
                for k in range(len(df2.iloc[j])):
                    percentage = df2.iloc[j][k] / sum_count * 100
                    df2.iloc[j][k] = percentage
            freq_list.append(df2.iloc[position-1])
        
    #create a graph
    df2 = pd.DataFrame(freq_list)
    df2 = df2.reset_index()
    df2 = df2.drop(df2.columns[0], axis=1)
    columns_list = df2.columns.values

    df2 = df2.values

    plt.style.use('seaborn-bright')
    length_f = df2.shape[1]

    bottom = 0

    fig, ax = plt.subplots()
    color = ['lightgray','red','blue','green','yellow','purple', 'cyan', 'magenta','pink', 'orange', 'gray', 'rosybrown', 'olive', 'deepskyblue', 'navy', 'yellowgreen', 'teal', 'moccasin', 'lime', 'salmon', 'aquamarine']
    for i in range(length_f):
        ax.bar(clade_names, df2[:, i], bottom=bottom, align="center", color=color[i], width=0.5)
        bottom += df2[:, i]

    ax.set_xlabel('clade')
    ax.set_ylabel('frequency')
    ax.legend(columns_list, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1.01,), title='amino acid', fontsize=5.5).get_title().set_fontsize(6.5)
    plt.title('{} {}'.format(protein, aa_change))

    plt.show()

    fig.savefig("/mnt/c/Users/naoki/Onedrive/デスクトップ/{}_{}.png".format(protein, aa_change), dpi=300, bbox_inches='tight')

input_protein = args.protein_name
input_aa_change = args.aa_change

make_graph_data(input_protein)
create_graph(input_protein, input_aa_change)

# %%
