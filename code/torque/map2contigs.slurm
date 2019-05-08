#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 150:00:00
#SBATCH --mem=170000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard


module load samtools

# Build bowtie database from contigs
/home/mj4pw/bin/bowtie2/bowtie2-build ${metagenome}.assembled/${metagenome}.final.contigs.1500.fasta ${metagenome}.assembled/${metagenome}_contigDB

# Align paired-end reads to contigs
/home/mj4pw/bin/bowtie2/bowtie2 -f --fr -x ${metagenome}.assembled/${metagenome}_contigDB \
                -1 ${metagenome}.R1.pool.cut.trim.fastq \
                -2 ${metagenome}.R1.pool.cut.trim.fastq \
                -S ${metagenome}.pair.reads2contigs.sam
samtools view -bS ${metagenome}.pair.reads2contigs.sam > ${metagenome}.pair.reads2contigs.bam
samtools sort ${metagenome}.pair.reads2contigs.bam ${metagenome}.pair.reads2contigs.sort
gzip ${metagenome}.pair.reads2contigs.sam
rm ${metagenome}.pair.reads2contigs.bam

# Align orphaned reads to contigs
/home/mj4pw/bin/bowtie2/bowtie2 -f -x ${metagenome}.assembled/${metagenome}_contigDB \
                -U ${metagenome}.pool.cut.trim.orphan.fastq \
                -S ${metagenome}.orphan.reads2contigs.sam
samtools view -bS ${metagenome}.orphan.reads2contigs.sam > ${metagenome}.orphan.reads2contigs.bam
samtools sort ${metagenome}.orphan.reads2contigs.bam ${metagenome}.orphan.reads2contigs.sort
gzip ${metagenome}.orphan.reads2contigs.sam
rm ${metagenome}.orphan.reads2contigs.bam

# Merge paired and orphan mappings
samtools merge ${metagenome}.reads2contigs.merge.bam ${metagenome}.pair.reads2contigs.sort.bam ${metagenome}.orphan.reads2contigs.sort.bam
samtools sort ${metagenome}.reads2contigs.merge.bam ${metagenome}.reads2contigs.merge.sort
rm ${metagenome}.reads2contigs.merge.bam ${metagenome}.pair.reads2contigs.sort.bam ${metagenome}.orphan.reads2contigs.sort.bam

# Remove PCR and optical duplicates
java -jar /home/mj4pw/bin/picard.jar MarkDuplicates \
	INPUT=${metagenome}.reads2contigs.merge.sort.bam \
	OUTPUT=${metagenome}.reads2contigs.merge.rmdup.bam \
	METRICS_FILE=${metagenome}.reads2contigs.merge.sort.rmdup.metrics \
	VALIDATION_STRINGENCY=LENIENT \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=6000 \
	REMOVE_DUPLICATES=TRUE
samtools sort ${metagenome}.reads2contigs.merge.rmdup.bam ${metagenome}.reads2contigs.merge.rmdup.sort
rm ${metagenome}.reads2contigs.merge.sort.bam ${metagenome}.reads2contigs.merge.rmdup.bam

