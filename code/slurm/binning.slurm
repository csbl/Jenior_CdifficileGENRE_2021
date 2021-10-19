#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 150:00:00
#SBATCH --mem=150000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard


# Generate depth file
/scratch/mj4pw/bin/metabat/jgi_summarize_bam_contig_depths --outputDepth ${metagenome}.depth.txt ${metagenome}.reads2contigs.merge.rmdup.sort.bam

# Bin contigs
gzip ${metagenome}.assembled/${metagenome}.final.contigs.1500.fasta
/scratch/mj4pw/bin/metabat/metabat2 -i ${metagenome}.assembled/${metagenome}.final.contigs.1500.fasta.gz \
									-a ${metagenome}.depth.txt \
									-o ${metagenome}_bins/bin

# Check bin quality
/scratch/mj4pw/bin/checkM/bin/checkm lineage_wf -t 8 -x fa ${metagenome}_bins ${metagenome}_bins/checkM

