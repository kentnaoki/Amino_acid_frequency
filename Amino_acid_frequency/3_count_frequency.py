from collections import Counter
import pandas as pd
import glob
import os

wd = "path to the folder to put files"
clade = ["G", "GR", "GH", "GV", "L", "O", "S", "V"]
protein = ["NS7a","NS7b","NS9b","NS9c","NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP11","NSP12","NSP13","NSP14","NSP15","NSP16","NS3", "NS6", "NS8", "N", "M", "Spike", "E"]

#create directory for each clade
try:
    os.mkdir(wd)
except:
    pass

try:
    for clade_name in clade:
        os.mkdir(wd + '/' + clade_name)
except:
    pass

gisaid_file = glob.glob("path to the directory which includes aligned sequence")
#open each file in the directory
for clade_file, clade_name in zip(gisaid_file, clade):
    protein_file = glob.glob('{}/*'.format(clade_file))
    for sep_file, protein_name in zip(protein_file, protein):

        #Count the number of amino acids in each position
        f = open(sep_file)
        len_seq = []
        seq = {}
        for line in f:
            if line.startswith('>'):
                namey = line[1:].rstrip()
                seq[namey] = ""
            else:
                seq[namey] += line.rstrip()
                len_seq.append(len(line.rstrip()))

        max_len = max(len_seq) #maximum length of sequences
        output = [''] * max_len
        for value in seq.values():
            for _id, c in enumerate(value):
                output[_id] += c

        #count amino acids
        count_aa = {}
        for posi in range(len(output)):
            count_aa[posi+1] = ''
            count_aa[posi+1] += str(Counter(output[posi]))

        df = pd.DataFrame(list(count_aa.items()), columns=['position', 'frequency'])
        df = df.set_index('position')
        #exclude unnecessary letters
        df['frequency'] = df['frequency'].str.replace('Counter', '')
        df['frequency'] = df['frequency'].str.replace('{', '')
        df['frequency'] = df['frequency'].str.replace('}', '')
        df['frequency'] = df['frequency'].str.replace('(', '')
        df['frequency'] = df['frequency'].str.replace(')', '')
        df['frequency'] = df['frequency'].str.replace("'", '')

        max_counter = []
        for num in range(1, df['frequency'].count()+1):
            max_counter.append(df['frequency'][num].count(','))
        aa_max = max(max_counter)+1

        tmp = df['frequency'].str.split(',|:', expand=True)
        i=0
        for tem_count in range(1, aa_max+1):
                df['a'+ str(tem_count)] = tmp[i]
                df['f'+ str(tem_count)] = tmp[i+1]
                i += 2
        df = df.drop('frequency', axis=1)

        df.to_csv("path to the csv file created")
        f.close()