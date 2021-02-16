# This scripts reads tab delimited files containing gene ids and numerical values (tpm) and create a matrix 
# It will also look for min number of replicate (no of replicates minus 1) has expression value higher than threshold 

import os
import sys
import pandas as pd
import csv

# sets expresion threshold
threshold = 1

# i = 0 is the script itself
first_accession = sys.argv[1]

name = os.path.splitext(os.path.basename(first_accession))[0]

matrix = {}

with open(first_accession) as f:
	ncs = f.read().splitlines()
	for x in ncs:	
		tra_id,value = x.split('\t')
		pair  = {}
		try:
			pair[name] = float(value)
			matrix[tra_id] = {}		
			matrix[tra_id][name] = float(value)
		except ValueError:
			continue

args = sys.argv

args_len = len(args)
# i = 1 is the first file
i = 2 

while i<args_len:
	name = os.path.splitext(os.path.basename(args[i]))[0]
	with open(args[i]) as f:
		ncs = f.read().splitlines()
		for x in ncs:
			tra_id,value = x.split('\t')
			pair  = {}
			try:
				pair[name] = float(value)
				if tra_id not in matrix:
					matrix[tra_id] = {}			
				matrix[tra_id][name]  = float(value)
			except ValueError:
				continue
	i = i+1		

df = pd.DataFrame(matrix).T

# adds 0 where data is missing
df.fillna(0, inplace=True)

# setting threshold for no of replicate with expression level one less than no of replicates 
replicate_threshold = args_len - 2

remove = []
for gene in matrix:
	count = 0
	for line in matrix[gene]:
		if matrix[gene][line] > threshold:
			count = count + 1
	
	if count < replicate_threshold:
		remove.append(gene)

# unique list of lncRNAs with less than threshold level expression			
remove = list(dict.fromkeys(remove))

# removing lowly expressed lncRNAs from matrix
for rem in remove:
	del matrix[rem]

df = pd.DataFrame(matrix).T

# replacing 0 with -
df.fillna(0, inplace=True)

# printing a tab-separated file
print(df.to_csv(sep='\t'))
