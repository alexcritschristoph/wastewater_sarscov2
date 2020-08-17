from Bio import SeqIO
from collections import defaultdict

genotypes = defaultdict(list)

import sys
correct = defaultdict(int)
for record in SeqIO.parse(sys.argv[1], 'fasta'):
	for pos in range(0, len(record.seq)):
		if record.seq[pos].upper() in ['A', 'C', 'G', 'T']:
			genotypes[pos].append(record.seq[pos].upper())

	if record.id == 'Wuhan/WIV04/2019':
		good = 0
		for pos in range(0, len(record.seq)):
			if record.seq[pos].upper() in ['A', 'C', 'G', 'T']:
				correct[pos] = good
				good += 1


import pandas as pd
for g in genotypes: #For each position
	if len(pd.Series(genotypes[g]).value_counts()) > 1: #if there was more than 1 nucleotide seen at this position
		passing_vars = 0
		sum = pd.Series(genotypes[g]).value_counts().sum()
		for i in pd.Series(genotypes[g]).value_counts():
			if float(i) >= int(sys.argv[2]): # if this allele was seen in at least N genomes
				passing_vars += 1
		if passing_vars > 1: 
			for var in pd.Series(genotypes[g]).value_counts().index:
				print(str(correct[g]) + "\t" + str(var) + "\t" + str(pd.Series(genotypes[g]).value_counts()[var] / pd.Series(genotypes[g]).value_counts().sum()))



