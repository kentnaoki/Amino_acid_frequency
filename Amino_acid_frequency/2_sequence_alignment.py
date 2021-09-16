import os
import glob

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
protein = ["E", "M", "N", "NS3", "NS6", "NS7a", "NS7b", "NS8", "NS9b", "NS9c", "NSP1", "NSP10", "NSP11", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP9", "Spike"]
try:
    os.mkdir(wd)
    for clade_name in clade:
        os.mkdir(wd + clade_name)
        for file_name in ["clustalo_in", "clustalo_out", "output"]:
            os.mkdir(wd + clade_name + "/" + file_name + "/")
except:
    pass

#align the sequence of FASTA file
gisaid_file = glob.glob("path to the directory having directories for each clade")
for clade_file, clade_name in zip(gisaid_file, clade):
    protein_file = glob.glob('{}/*'.format(clade_file))
    for in_file, out_file in zip(protein_file, protein):
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
            newf = open(wd + clade_name + "/output/" + clade_name + "_" + out_file + ".aln", "w")
            for k, v in aln_in.items():
                for k2, v2 in sequence_dict.items():
                    if k in v2:
                        for x in v2:
                            newf.write(x + "\n")
                            newf.write(v + "\n")
            newf.close()

        aln_in = aln2dict(wd + clade_name + "/clustalo_out/" + clade_name + "_" + out_file + "_clustalo" + ".fasta")

        repopulate(aln_in, sequence_dict)