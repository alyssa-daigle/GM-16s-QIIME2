#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=512G
#SBATCH --job-name="cut-dada"
#SBATCH --output="job_%j.log"
#SBATCH --mail-type=END
#SBATCH --mail-user=and1038@usnh.edu
#SBATCH --exclude=node\[117-118\]

start=/mnt/gpfs01/home/obrien/and1038/GMproject/qiimeoutputs/

## Load the necessary environment
module purge
module load anaconda/colsa
conda activate qiime2-2024.10.26


echo "starting data import at" $(date)

qiime tools import \
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path /mnt/gpfs01/home/obrien/and1038/GMproject/manifest.csv \
 --output-path ${start}GM-paired-end.qza \
 --input-format PairedEndFastqManifestPhred33

echo "done import at" $(date)


echo "starting demux summary at" $(date)

qiime demux summarize \
  --i-data ${start}GM-paired-end.qza \
  --o-visualization $start}GM-paired-end.qzv

echo "done demux summary at" $(date)

echo "starting cutadapt at" $(date)

#replace with primers used in your specific analysis
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences ${start}GM-paired-end.qza \
 --p-front-f GTGYCAGCMGCCGCGGTAA \ 
 --p-front-r CCGYCAATTYMTTTRAGTTT \
 --o-trimmed-sequences ${start}GM-trimmed.qza

echo "done cutadapt at"  $(date)


echo "starting output (trimmed) visualization at" $(date)

qiime demux summarize \
 --i-data ${start}GM-trimmed.qza \
 --o-visualization ${start}GM-trimmed.qzv

echo "done trimmed output visualization at" $(date)


echo "starting denoise at" $(date)

#trunc length here is based on the demux results (GM-trimmed.qzv) which tells us where quality beings to drop off for fwd and rev reads
#quality score stays high through 240 bp for fwd reads and drops off around 200 bp for rev reads
#want at least 20 bp of overlap here
#since the V4-V5 regions of 16s rRNA gene is about 400 bp in length, the total read length after truncation is 240 + 200 bp = 440
#400bp 16s rRNA - (240 fwd + 200 rev) = 40 bp overlap, which is good

qiime dada2 denoise-paired \
 --i-demultiplexed-seqs ${start}GM-trimmed.qza \
 --p-trunc-len-f 240 \
 --p-trunc-len-r 200 \
 --p-pooling-method 'pseudo' \
 --p-n-threads 24 \
 --o-representative-sequences ${start}dada-repseqs.qzv \
 --o-table ${start}dada-table.qza \
 --o-denoising-stats ${start}dada-stats.qza

echo "done denoise at" $(date)


echo "starting dada stats output visualization at" $(date)

qiime metadata tabulate \
 --m-input-file ${start}dada-stats.qza \
 --o-visualization ${start}dada-stats.qzv

 echo "done dada stats output vis at" $(date)


echo "starting dada (feature) table vis at" $(date)

qiime feature-table summarize \
 --i-table ${start}dada-table.qza \
 --o-visualization ${start}dada-table.qzv \
 --m-sample-metadata-file /mnt/gpfs01/home/obrien/and1038/GMproject/metadata.tsv

echo "done dada (feature) table vis at" $(date)


echo "starting dada rep seqs vis at" $(date)

qiime feature-table tabulate-seqs \
 --i-data ${start}dada-repseqs.qza \
 --o-visualization ${start}dada-repseqs.qzv

echo "done dada rep seqs vis at" $(date)


echo "starting export of feature table at" $(date)

#outputs feature table to use later in R
qiime tools export \
 --input-path ${start}dada-table.qza \
 --output-path ${start}feature-table

 #biom convert -i ${start}feature-table/feature-table.biom -o ${start}feature-table/feature-table.tsv --to-tsv

 echo "done export of feature table at" $(date)
