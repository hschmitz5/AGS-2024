#!/bin/bash
#SBATCH --job-name="import"
#SBATCH --output=import.out
#SBATCH --error=import.err
#SBATCH -A p31629
#SBATCH -p normal
#SBATCH -t 00:15:00
#SBATCH -N 1
#SBATCH --mem=1G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hannahschmitz2026@u.northwestern.edu

# activate QIIME2 
module purge all
module load qiime2/2025.4-amplicon

# import reads into qiime-comatible format
qiime tools import \
--input-path /projects/b1052/hannah/data_2930/qiime_io/manifest.txt \
--output-path /projects/b1052/hannah/data_2930/qiime_io/reads.qza \
--input-format PairedEndFastqManifestPhred33V2 \
--type SampleData[PairedEndSequencesWithQuality]

qiime demux summarize \
--i-data /projects/b1052/hannah/data_2930/qiime_io/reads.qza \
--o-visualization /projects/b1052/hannah/data_2930/qiime_io/readquality_raw.qzv

## using a space and backslash allows you to insert a linebreak without disrupting the function
## you can also have the command written as one line but this is harder to read
