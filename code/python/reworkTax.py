#!/usr/bin/python
'''USAGE: reworkTax.py taxonomy_file
This script reformats mothur taxonomy files
'''
import sys

outFile = str(sys.argv[1]).replace('taxonomy','') + 'format.taxonomy'
outFile = open(outFile, 'w')
outFile.write('otu\tphylum\tclass\torder\tfamily\tgenus\n')

with open(sys.argv[1], 'r') as taxonomy:

	for line in taxonomy:

		line = line.split()
		if line[0] == 'OTU':
			continue

		curr_otu = line[0]
		curr_tax = line[2].split(';')
		curr_phylum = curr_tax[1]
		curr_class = curr_tax[2]
		curr_order = curr_tax[3]
		curr_family = curr_tax[4]
		curr_genus = curr_tax[5]

		entry = curr_otu + '\t' + curr_phylum + '\t' + curr_class + '\t' + curr_order + '\t' + curr_family + '\t' + curr_genus + '\n'
		outFile.write(entry)

outFile.close()

