#!/bin/bash
#SBATCH --job-name="taxonomy"
#SBATCH --output=taxonomy.out
#SBATCH --error=taxonomy.err
#SBATCH -A p31629
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hannahschmitz2026@u.northwestern.edu

module purge all
module load qiime2/2025.4-amplicon

# classify sequences from previous step
## keep the classifier filepath the same
qiime feature-classifier classify-sklearn \
--i-classifier /projects/b1052/hannah/resources/midas_5.3_classifier.qza \
--i-reads /projects/b1052/hannah/data_2930/qiime_io/rep_seqs_dada2.qza \
--o-classification /projects/b1052/hannah/data_2930/qiime_io/taxonomy.qza

qiime metadata tabulate \
--m-input-file /projects/b1052/hannah/data_2930/qiime_io/taxonomy.qza \
--o-visualization /projects/b1052/hannah/data_2930/qiime_io/taxonomy.qzv
