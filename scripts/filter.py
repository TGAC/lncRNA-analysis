# This scripts reads two tab delimited files containing gene ids and transcript ids and looks for non-matching transcripts in of file 2

import sys

nc_transcript = sys.argv[1] #"/hpc-home/thankia/myProject/GFF/nc_transcript/ARI_D.tsv2"
transcript = sys.argv[2] # "/hpc-home/thankia/myProject/GFF/transcript/ARI_D.tsv"

with open(nc_transcript) as f:
	ncs = f.read().splitlines() 

with open(transcript) as f:
	tra = f.read().splitlines() 

#Reverses array to store transcript id - gene id pairs
all_gene = dict(x.split('\t')[::-1] for x in tra)

#stores transcript ids
ncs_tra=[]
for x in ncs:
	ncs_tra.append(x.split('\t')[1])

all_tra=[]
for x in tra:
	all_tra.append(x.split('\t')[1])

#find the missing transcripts in ncs_tra
set_difference = set(all_tra) - set(ncs_tra)

#converts to list
list_difference = list(set_difference)


diff_gene = []

#stores gene ids of missing genes
for diff in list_difference:
	diff_gene.append(all_gene[diff])

#print(list(set(diff_gene)))
# using set removes duplicates
print("\n".join(set(diff_gene)))

