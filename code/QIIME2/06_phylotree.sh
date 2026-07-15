#!/bin/bash
#SBATCH --job-name="tree"
#SBATCH --output=tree.out
#SBATCH --error=tree.err
#SBATCH -A p31629
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hannahschmitz2026@u.northwestern.edu

module purge all
module load qiime2/2025.4-amplicon

# make phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences /projects/b1052/hannah/data_2930/qiime_io/rep_seqs_dada2.qza \
--o-alignment /projects/b1052/hannah/data_2930/qiime_io/rep_seqs_dada2_aligned.qza \
--o-masked-alignment /projects/b1052/hannah/data_2930/qiime_io/rep_seqs_dada2_aligned_masked.qza \
--o-tree /projects/b1052/hannah/data_2930/qiime_io/unrooted_tree.qza \
--o-rooted-tree /projects/b1052/hannah/data_2930/qiime_io/rooted_tree.qza
