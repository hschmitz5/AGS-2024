#!/bin/bash
#SBATCH --job-name="denoise"
#SBATCH --output=denoise.out
#SBATCH --error=denoise.err
#SBATCH -A p31629
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hannahschmitz2026@u.northwestern.edu

module purge all
module load qiime2/2025.4-amplicon


# run dada2 to identify ASVs
qiime dada2 denoise-paired --verbose \
--p-trunc-len-f 240 --p-trunc-len-r 224 \
--i-demultiplexed-seqs /projects/b1052/hannah/data_2930/qiime_io/reads_trimmed.qza \
--o-representative-sequences /projects/b1052/hannah/data_2930/qiime_io/rep_seqs_dada2.qza \
--o-table /projects/b1052/hannah/data_2930/qiime_io/table_dada2.qza \
--o-denoising-stats /projects/b1052/hannah/data_2930/qiime_io/stats_dada2.qza

## p trunc len f and r should be determined based on your data quality

# make visualization files
qiime metadata tabulate \
--m-input-file /projects/b1052/hannah/data_2930/qiime_io/stats_dada2.qza \
--o-visualization /projects/b1052/hannah/data_2930/qiime_io/stats_dada2.qzv

qiime feature-table tabulate-seqs \
--i-data /projects/b1052/hannah/data_2930/qiime_io/rep_seqs_dada2.qza \
--o-visualization /projects/b1052/hannah/data_2930/qiime_io/rep_seqs_dada2.qzv
