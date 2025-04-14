#!/bin/bash

cd /mnt/gpfs01/home/obrien/and1038/GMproject/taxid_filter

out=/mnt/gpfs01/home/obrien/and1038/GMproject/taxid_filter/

module purge
module load anaconda/colsa
conda activate ncbi_datasets


##make script executable if not done already

#chmod +x taxid_filter.sh

##command that runs script

# nohup bash /mnt/gpfs01/home/obrien/and1038/GMproject/taxid_filter/taxid_filter.sh \
# > /mnt/gpfs01/home/obrien/and1038/GMproject/taxid_filter/taxid_filter.out \
# 2> /mnt/gpfs01/home/obrien/and1038/GMproject/taxid_filter/taxid_filter.err &

##check status 

#ps aux | grep taxid_filter.sh



echo "starting Streptophyta"

datasets summary taxonomy taxon 35493 --children --as-json-lines | \
dataformat tsv taxonomy --template tax-summary > ${out}ncbi_taxid_Streptophyta.txt

echo "done Streptophyta"

echo "starting Chordata"

datasets summary taxonomy taxon 7711 --children --as-json-lines | \
dataformat tsv taxonomy --template tax-summary > ${out}ncbi_taxid_Chordata.txt

echo "done Chordata"

