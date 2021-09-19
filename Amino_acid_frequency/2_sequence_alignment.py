import os
import glob
import re

def fasta2dict(seq):
	#convert seqeuence to fasta
	clus = {}
	with open(seq) as files:
		for line in files:
			if line.startswith(">"):
				namey = line.rstrip()
				clus[namey] = []
			else:
				clus[namey].append(line.rstrip())
	for k, v in clus.items():
		clus[k] = "".join(v)
	return clus

def human_only(clus):
	clus2 = {}
	for k, v in clus.items():
		clus2[k] = v
	return clus2

def filter_seqs2(clus):
	# remove all non-unique sequences 
	sequence_dict = {}
	for seq_id, sequence in clus.items():
		if not sequence in sequence_dict:
			sequence_dict[sequence] = []
		sequence_dict[sequence].append(seq_id[:])
	for k in list(sequence_dict.keys()):
		if len(sequence_dict[k]) == 1:
			sequence_dict.pop(k)
	return sequence_dict

def round_down(num, divisor):
    return num - (num%divisor)

def make_clastalo_in(sequence_dict, namey):
	# make clastalo input 
	newf = open(namey , "w")
	seq_size = []
	for k, v in sequence_dict.items():
		seq_size.append(float(len(k)))

	avg_size = sum(seq_size)/len(seq_size)
	avg_size =  round_down(avg_size,10)


	for k, v in sequence_dict.items():
		if not "X" in k and len(k) >= avg_size:
			newf.write(v[0] + "\n")
			newf.write(k + "\n")
	newf.close()

def clustalo(alnfile, outfile):
	os.system("/usr/bin/clustalo -i %s -o %s "%(alnfile, outfile))

#create directory to put files
wd = "path to the directory"
clade = ["L", "S", "V", "G", "GH", "GR", "GV", "O"]

try:
    os.mkdir(wd)
    for clade_name in clade:
        os.mkdir(wd + clade_name)
        for file_name in ["clustalo_in", "clustalo_out", "output"]:
            os.mkdir(wd + clade_name + "/" + file_name + "/")
except:
    pass

#align the sequence of FASTA file
gisaid_file = glob.glob("path to the FASTA files")
for in_file in gisaid_file:
    
    clade_name = re.split('[/_.]', in_file)[-3]
    protein = re.split('[/_.]', in_file)[-2]
    clus = fasta2dict(in_file)
    clus = human_only(clus)
    sequence_dict = filter_seqs2(clus)
    make_clastalo_in(sequence_dict, wd + "path to FASTA file")

    clustalo(wd + "path to FASTA file", wd + "path to FASTA file")

    def aln2dict(namey):
        # Convert alignment to dictionary 
        aln_in = {}
        with open(namey) as aln:
            for line in aln:
                if line.startswith(">"):
                    namey = line.rstrip()
                    aln_in[namey] = []
                else:
                    aln_in[namey].append(line.rstrip())
        aln2 = {}
        for k, v in aln_in.items():
            aln2[k[:]] = "".join(v)
        return aln2

    def repopulate(aln_in, sequence_dict):
        # repopulate alignment with non-unique sequences 
        newf = open(wd + clade_name + "/output/" + clade_name + "_" + protein + ".aln", "w")
        for k, v in aln_in.items():
            for k2, v2 in sequence_dict.items():
                if k in v2:
                    for x in v2:
                        newf.write(x + "\n")
                        newf.write(v + "\n")
        newf.close()

    aln_in = aln2dict(wd + clade_name + "/clustalo_out/" + clade_name + "_" + protein + "_clustalo" + ".fasta")

    repopulate(aln_in, sequence_dict)