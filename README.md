This repository contains the complete pipeline for the 16s analysis of my master's research, the Green Manures project. Scripts are labeled numerically in the order they are used in the pipeline. Both [QIIME2](https://docs.qiime2.org/2024.10/tutorials/qiime2-for-experienced-microbiome-researchers/) and R are used for data processing. 

* **qiime_inputs** - Directory containing necessary inputs for qiime scripts. metadata.tsv contains the metadata for my 263 samples. manifest.csv contains the absolute filepaths of my raw reads. 
* **R_inputs** - Directory containing necessary inputs for R scripts. Feature table not included due to large file size.

* All sequence data goes through scripts 01-06, and Expt 1 and Expt 2 analyses are then separate:
    * * **01_cutadapt_dada2_GM16s.sh** - Slurm script for data importing, cutadapt, trimming, and denoising in qiime.
    * **02_CyanoSeq_GM16s.sh** - Slurm script for training a classifier and assigning taxonomy. I am using the [CyanoSeq + SILVA combined classifier](https://zenodo.org/records/13910424) as I am looking for higher resolution when it comes to cyanobacterial reads.
  
**NOTE**: I chose to do all my filtering of unwanted reads/samples in R for more control, but it can be done in QIIME as well.
  
    * **03_taxize.R**  - R script using the [taxize](https://github.com/ropensci/taxize) package to take the assigned taxonomy from script 02 and assign taxa IDs to families using the NCBI database. This step is necessary for later filtering out unwanted reads (plant or human-associated). 
    * **04_taxid_StrepChord_filter.sh** - Slurm script for obtaining taxa IDs from NCBI for taxa that I will eventually be removing from my data (Chordata and Streptophyta).
    * **05_data_formatting.R** - R script that formats tables.
    * **06_filter_taxaIDs.R** - R script that filters data based on taxa IDs determined in the previous scripts. Blanks are also removed at this step.

**NOTE**: The command `qiime tools export` exports the feature table as a .biom file in script 01, so I used `biom convert -i feature-table.biom -o feature-table.tsv --to-tsv` to convert the table to a tsv to be used in R.

