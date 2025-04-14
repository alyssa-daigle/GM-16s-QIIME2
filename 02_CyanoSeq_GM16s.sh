#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=512G
#SBATCH --job-name="taxacoll"
#SBATCH --output="job_%j.log"
#SBATCH --mail-type=END
#SBATCH --mail-user=and1038@usnh.edu
#SBATCH --exclude=node\[117-118\]

start=/mnt/gpfs01/home/obrien/and1038/GMproject/qiimeoutputs/

## Load the necessary environment
module purge
module load anaconda/colsa
conda activate qiime2-2024.10.26

# using SILVA + CyanoSeq database
echo "starting CyanoSeq data import at" $(date)

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /mnt/gpfs01/home/obrien/and1038/GMproject/training-feature-classifiers/CyanoSeqv1.3_SILVA138.2_sequences.fasta \
  --output-path /mnt/gpfs01/home/obrien/and1038/GMproject/training-feature-classifiers/CyanoSeq-ref-seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /mnt/gpfs01/home/obrien/and1038/GMproject/training-feature-classifiers/CyanoSeqv1.3_SILVA138.2_taxonomy.tsv \
  --output-path /mnt/gpfs01/home/obrien/and1038/GMproject/training-feature-classifiers/CyanoSeq-ref-taxonomy.qza

echo "done CyanoSeq data import at" $(date)


echo  "starting extraction of region-specific reads at" $(date)

qiime feature-classifier extract-reads \
 --i-sequences /mnt/gpfs01/home/obrien/and1038/GMproject/training-feature-classifiers/CyanoSeq-ref-seqs.qza \
 --p-f-primer GTGYCAGCMGCCGCGGTAA \ 
 --p-r-primer CCGYCAATTYMTTTRAGTTT \
 --p-trunc-len 490 \
 --p-min-length 231 \
 --o-reads /mnt/gpfs01/home/obrien/and1038/GMproject/training-feature-classifiers/CyanoSeq-ref-seqs-trimmed.qza

echo  "done read extraction at" $(date)


echo "starting CyanoSeq classifier training at" $(date)

qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads /mnt/gpfs01/home/obrien/and1038/GMproject/training-feature-classifiers/CyanoSeq-ref-seqs-trimmed.qza \
 --i-reference-taxonomy /mnt/gpfs01/home/obrien/and1038/GMproject/training-feature-classifiers/CyanoSeq-ref-taxonomy.qza \
 --o-classifier /mnt/gpfs01/home/obrien/and1038/GMproject/training-feature-classifiers/CyanoSeq-trained-classifier.qza

echo "done classifier training at" $(date)


echo "starting taxonomy assignment at" $(date)

qiime feature-classifier classify-sklearn \
 --i-classifier /mnt/gpfs01/home/obrien/and1038/GMproject/training-feature-classifiers/CyanoSeq-trained-classifier.qza \
 --i-reads ${start}dada-repseqs.qza \
 --p-confidence 0.7 \
 --o-classification ${start}taxonomy.qza \
 --p-n-jobs 0

echo "done taxonomy assignment at" $(date)


echo "starting taxonomy assign output at" $(date)

qiime metadata tabulate \
 --m-input-file ${start}taxonomy.qza \
 --o-visualization ${start}taxonomy.qzv

echo "done taxonomy assign ouput at" $(date)


echo "starting export of taxonomy as tsv at" $(date)

qiime tools export \
  --input-path ${start}taxonomy.qza \
  --output-path ${start}taxonomy

echo "done export of taxonomy as tsv at" $(date)


echo "starting barplot at" $(date)

qiime taxa barplot \
 --i-table ${start}noplants-feature-table.qza \
 --i-taxonomy ${start}taxonomy.qza \
 --m-metadata-file /mnt/gpfs01/home/obrien/and1038/GMproject/metadata.tsv \
 --o-visualization ${start}prefilter_barplot.qzv

echo "done barplot at" $(date)




