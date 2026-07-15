#!/bin/bash
#SBATCH --job-name="classifier"
#SBATCH --output=classifier.out
#SBATCH --error=classifier.err
#SBATCH -A p31629
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=20G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hannahschmitz2026@u.northwestern.edu

module purge all
module load qiime2/2025.4-amplicon

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /projects/b1052/hannah/resources/midas_5.3.fa \
  --output-path /projects/b1052/hannah/resources/midas_seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /projects/b1052/hannah/resources/midas_5.3.txt \
  --output-path /projects/b1052/hannah/resources/midas_tax.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /projects/b1052/hannah/resources/midas_seqs.qza \
  --i-reference-taxonomy /projects/b1052/hannah/resources/midas_tax.qza \
  --o-classifier /projects/b1052/hannah/resources/midas_5.3_classifier.qza
