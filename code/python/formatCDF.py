#!/usr/bin/python
import sys
import numpy

with open(sys.argv[1], 'r') as inFile:

	outFile = open(sys.argv[2], 'w')
	header = inFile.readline()
	header = header.split()
	header = ','.join(header) + '\n'
	outFile.write(header)

	for line in inFile:
		line = line.split()
		transcript = []
		for x in line[1:]:
			x = float(x)
			#x = numpy.log2(x + 1.0)
			x = str(x)
			transcript.append(x)

		transcript = ','.join(transcript)
		genes = line[0].split(',')

		for gene in genes:
			entry = gene + ',' + transcript + '\n'
			outFile.write(entry)

outFile.close()
