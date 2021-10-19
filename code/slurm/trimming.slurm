#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 150:00:00
#SBATCH --mem=90000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard


# Assign correct Nextera XT primer pair that was used (5' + 3' reverse complement)
if [ $sample = 'noAbx_mock' ]
then
	adapt_F=CAAGCAGAAGACGGCATACGAGATTAAGGCGAGTCTCGTGGGCTCGG
	adapt_R=GACGCTGCCGACGAGCGATCTAGTGTAGATCTCGGTGGTCGCCGTATCATT
elif [ $sample = 'cefoperazone_mock' ]
then
	adapt_F=CAAGCAGAAGACGGCATACGAGATAGGCAGAAGTCTCGTGGGCTCGG
	adapt_R=GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT
elif [ $sample = 'streptomycin_mock' ]
then
	adapt_F=CAAGCAGAAGACGGCATACGAGATTCCTGAGCGTCTCGTGGGCTCGG
	adapt_R=GACGCTGCCGACGATCTACTCTGTGTAGATCTCGGTGGTCGCCGTATCATT
elif [ $sample = 'clindamycin_mock' ]
then
	adapt_F=CAAGCAGAAGACGGCATACGAGATCGTACTAGGTCTCGTGGGCTCGG
	adapt_R=AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC
fi

# Cut residual Illumina adapters
/home/mj4pw/.local/bin/cutadapt --error-rate=0.1 --overlap=10 \
	-a ${adapt_F} \
	-A ${adapt_R} \
	-o ${sample}.R1.pool.cut.fastq \
	-p ${sample}.R2.pool.cut.fastq \
	${sample}.R1.pool.fastq \
	${sample}.R2.pool.fastq


# Quality trim the samples based on quality scores
/home/mj4pw/bin/sickle/sickle pe -t sanger -q 25 -l 21 \
	-f ${sample}.R1.pool.cut.fastq \
	-r ${sample}.R2.pool.cut.fastq \
	-o ${sample}.R1.pool.cut.trim.fastq \
	-p ${sample}.R2.pool.cut.trim.fastq \
	-s ${sample}.pool.cut.trim.orphan.fastq 


# Record stats
python /scratch/mj4pw/bin/analysis-scripts/readStats.py ${sample}.R1.pool.cut.trim.fastq > ${sample}.trimming.txt
python /scratch/mj4pw/bin/analysis-scripts/readStats.py ${sample}.R2.pool.cut.trim.fastq >> ${sample}.trimming.txt
python /scratch/mj4pw/bin/analysis-scripts/readStats.py ${sample}.cut.trim.orphan.fastq >> ${sample}.trimming.txt

