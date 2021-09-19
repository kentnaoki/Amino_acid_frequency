from collections import Counter
import pandas as pd
import glob
import os
import re

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

#function to count the amino acids in each position
def make_freq_file():
    gisaid_file = glob.glob("path to the directory which includes aligned sequence")

    #Count the number of amino acids in each position
    for sep_file in gisaid_file:
        protein_name = re.split('[/_.]', sep_file)[-2]
        clade_name = re.split('[/_.]', sep_file)[-3]
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

        df.to_csv("path to the directory where you want to put the file/{}_{}.csv".format(clade_name, protein_name))
        f.close()

#function to check the positions that the most frequent amino acids are different
def check_most_freq():
    freq_files = glob.glob("The directory where you put the csv files")
    for freq_file in freq_files:
        protein_name = re.split('[/_.]', freq_file)[-2]
        df = pd.read_csv('{}'.format(freq_file))
    
        df2 = df[~(df['L'] == df['S']) | ~(df['S'] == df['V']) | ~(df['V'] == df['O']) | ~(df['O'] == df['G']) | ~(df['G'] == df['GH']) | ~(df['GH'] == df['GR']) | ~(df['GR'] == df['GV'])]
        if len(df2.index) == 0:
            pass
        else:
            df2.to_csv('path to the directory you want to put the files/{}.csv'.format(protein_name), index=False)

make_freq_file()
check_most_freq()