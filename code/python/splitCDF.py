#!/usr/bin/python
'''USAGE: splitCDF.py
This script splits cdf in vivo transcription data and creates python dictionarie for each condition
'''
import sys
import pickle

strep_dict = {}
cef_dict = {}
clinda_dict = {}
gf_dict = {}
gene_set = set()

with open(sys.argv[1], 'r') as transcript:

	for line in transcript:
		line = line.split()
		if line[0] == 'KO':
			continue

		cef_curr = int(line[1])
		clinda_curr = int(line[2])
		strep_curr = int(line[3])
		gf_curr = int(line[4])
		genes_curr = line[6].split(',')
		
		for index in genes_curr:
			if not index in gene_set:
				gene_set.add(index)
				cef_dict[index] = cef_curr
				clinda_dict[index] = clinda_curr
				strep_dict[index] = strep_curr
				gf_dict[index] = gf_curr
			else:
				cef_dict[index] = cef_dict[index] + cef_curr
				clinda_dict[index] = clinda_dict[index] + clinda_curr
				strep_dict[index] = strep_dict[index] + strep_curr
				gf_dict[index] = gf_dict[index] + gf_curr


pickle.dump(cef_dict, open('cef_dict.pickle', 'wb'))
pickle.dump(clinda_dict, open('clinda_dict.pickle', 'wb'))
pickle.dump(strep_dict, open('strep_dict.pickle', 'wb'))
pickle.dump(gf_dict, open('gf_dict.pickle', 'wb'))



