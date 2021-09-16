#%%
import pandas as pd
import matplotlib.pyplot as plt
import glob

#align data
gisaid_file = glob.glob("path to the file having the data for creating graph")

protein_names = ['N', 'NS9c', 'NSP12', 'NSP2']
clade_names = ['G', 'GH', 'GR', 'GV', 'L', 'O', 'S', 'V']

N_203_list = []
N_220_list = []
NS9c_54_list = []
NS9c_71_list = []
NSP12_list = []
NSP2_list = []
nested_list = []
for protein_folder, protein_name in zip(gisaid_file, protein_names):
    clade_folder = glob.glob('{}/*'.format(protein_folder))
    for folder, clade_name in zip(clade_folder, clade_names):
        graph_file = glob.glob('{}/*'.format(folder))
        graph_file = graph_file[0]
        
        if graph_file[-8:-4] == 'NSP2':
            df = pd.read_csv(graph_file)
            df2 = df.drop(df.columns[0], axis=1)
            df2 = df2.fillna(0).astype('int64', errors='ignore')
            df2 = df2.astype(float)
#calculate the ratio of amino acids in NSP2
            sum_count = 0
            for m in range(len(df2.iloc[0])):
                sum_count += df2.iloc[0][m]
            for j in range(len(df2.index)):
                for k in range(len(df2.iloc[j])):
                    percentage = df2.iloc[j][k] / sum_count * 100
                    df2.iloc[j][k] = percentage
            NSP2_list.append(df2.iloc[84])

        elif protein_name == 'N':
            df = pd.read_csv(graph_file)
            df2 = df.drop(df.columns[0], axis=1)
            df2 = df2.fillna(0).astype('int64', errors='ignore')
            df2 = df2.astype(float)         
#calculate the ratio of amino acids in N protein
            sum_count = 0
            for m in range(len(df2.iloc[0])):
                sum_count += df2.iloc[0][m]
            for j in range(len(df2.index)):
                for k in range(len(df2.iloc[j])):
                    percentage = df2.iloc[j][k] / sum_count * 100
                    df2.iloc[j][k] = percentage
            N_203_list.append(df2.iloc[202])
            N_220_list.append(df2.iloc[219])

        elif protein_name == 'NS9c':
            df = pd.read_csv(graph_file)
            df2 = df.drop(df.columns[0], axis=1)
            df2 = df2.fillna(0).astype('int64', errors='ignore')
            df2 = df2.astype(float)       
#calculate the ratio of amino acids in NS9c
            sum_count = 0
            for m in range(len(df2.iloc[0])):
                sum_count += df2.iloc[0][m]
            for j in range(len(df2.index)):
                for k in range(len(df2.iloc[j])):
                    percentage = df2.iloc[j][k] / sum_count * 100
                    df2.iloc[j][k] = percentage
            NS9c_54_list.append(df2.iloc[53])
            NS9c_71_list.append(df2.iloc[70])

        elif protein_name == 'NSP12':
            df = pd.read_csv(graph_file)
            df2 = df.drop(df.columns[0], axis=1)
            df2 = df2.fillna(0).astype('int64', errors='ignore')
            df2 = df2.astype(float)       
#calculate the ratio of amino acids in NSP12
            sum_count = 0
            for m in range(len(df2.iloc[0])):
                sum_count += df2.iloc[0][m]
            for j in range(len(df2.index)):
                for k in range(len(df2.iloc[j])):
                    percentage = df2.iloc[j][k] / sum_count * 100
                    df2.iloc[j][k] = percentage
            NSP12_list.append(df2.iloc[322])

#append to nested list
nested_list.append(N_203_list)
nested_list.append(N_220_list)
nested_list.append(NS9c_54_list)
nested_list.append(NS9c_71_list)
nested_list.append(NSP2_list)
nested_list.append(NSP12_list)

proteins = ['N', 'N', 'NS9c', 'NS9c', 'NSP2', 'NSP12']
positions = ['R203K','A220V', 'G54N', 'L71F', 'T85I', 'P323L']

#create a graph
for protein_list, protein, position in zip(nested_list, proteins, positions):

    df2 = pd.DataFrame(protein_list)
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
    plt.title('{} {}'.format(protein, position))

    plt.show()

    fig.savefig("path and the name of the figure", dpi=300, bbox_inches='tight')
# %%
