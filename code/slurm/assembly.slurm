#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 150:00:00
#SBATCH --mem=150000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard


# Run Megahit to assemble contigs
python /home/mj4pw/bin/megahit/megahit -t 16 -m 45000000000 -l 250 \
	-1 ${metagenome}.R1.pool.cut.fastq \
	-2 ${metagenome}.R2.pool.cut.fastq \
	-r ${metagenome}.pool.cut.trim.orphan.fastq \
	--k-min 87 --k-max 127 --k-step 10 -o ${metagenome}.assembled

cd ${metagenome}.assembled

# Remove contigs shorter than read length
python /scratch/mj4pw/bin/analysis-scripts/rmShortSeqs.py final.contigs.fa ${metagenome}.final.contigs.1500.fa 1500

# Return some stats on the final contigs
python /scratch/mj4pw/bin/analysis-scripts/contigStats.py ${metagenome}.final.contigs.1500.fa > ${metagenome}.assembly.txt




