#!/usr/bin/python
'''USAGE: codonUsage.py genes.fasta codon_useage_table
Calculates codon utilization for each gene in a fasta file
'''
import sys
import pickle

codons = ['UUU', 'UCU', 'UAU', 'UGU', 'UUC', 'UCC', 'UAC', 'UGC', 'UUA', 'UCA', 'UAA', 'UGA', 'UUG', 'UCG', 'UAG', 'UGG', 'CUU', 'CCU', 'CAU', 'CGU', 'CUC', 'CCC', 'CAC', 'CGC', 'CUA', 'CCA', 'CAA', 'CGA', 'CUG', 'CCG', 'CAG', 'CGG', 'AUU', 'ACU', 'AAU', 'AGU', 'AUC', 'ACC', 'AAC', 'AGC', 'AUA', 'ACA', 'AAA', 'AGA', 'AUG', 'ACG', 'AAG', 'AGG', 'GUU', 'GCU', 'GAU', 'GGU', 'GUC', 'GCC', 'GAC', 'GGC', 'GUA', 'GCA', 'GAA', 'GGA', 'GUG', 'GCG', 'GAG', 'GGG']
outFile = open(sys.argv[2],'w')
entry = ['gene'] + codons
entry = '\t'.join(entry) + '\n'
outFile.write(entry)


with open(sys.argv[1],'r') as gene_seqs:

	gene = gene_seqs.readline()
	gene = gene.split()[0]
	gene = gene.replace('>cdf:','')
	seq = ''

	for line in gene_seqs:

		if line[0] == '>':
			seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
			codon_dict = {}
			for codon in seq:
				if len(codon) < 3:
					continue
				try:
					codon_dict[codon] = codon_dict[codon] + 1
				except KeyError:
					codon_dict[codon] = 1
			
			entry = [gene]
			for codon in codons:
				try:
					entry.append(str(codon_dict[codon]))
				except KeyError:
					entry.append('0')

			entry = '\t'.join(entry) + '\n'
			outFile.write(entry)

			gene = line.split()[0]
			gene = gene.replace('>cdf:','')
			seq = ''
			continue

		else:
			line = line.strip()
			line = line.upper()
			line = line.replace('T','U')
			seq = seq + line


	seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
	codon_dict = {}
	for codon in seq:
		if len(codon) < 3:
			continue
		try:
			codon_dict[codon] = codon_dict[codon] + 1
		except KeyError:
			codon_dict[codon] = 1
	
	entry = [gene]
	for codon in codons:
		try:
			entry.append(str(codon_dict[codon]))
		except KeyError:
			entry.append('0')

	entry = '\t'.join(entry) + '\n'
	outFile.write(entry)


outFile.close()
