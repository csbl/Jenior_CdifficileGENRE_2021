
import sys

outFasta = str(sys.argv[1]).rstrip('fastn') + 'format.fasta'
outFasta = open(outFasta, 'w')

with open(sys.argv[1], 'r') as inFasta:
	entry = inFasta.readline()
	entry = entry.split('|')[1] + '|'
	outFasta.write(entry)
	seq = ''

	for line in inFasta:
		if line[0] == '>':
			outFasta.write(seq + '\n')
			seq = ''
			entry = line.split('|')[1] + '|'
			outFasta.write(entry)
		else:
			seq += line.strip().upper()

outFasta.write(seq + '\n')
outFasta.close()
