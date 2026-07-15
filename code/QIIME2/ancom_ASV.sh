#!/bin/bash
#SBATCH --job-name="ancom_ASV"
#SBATCH --output=ancom_ASV-%j.out
#SBATCH --error=ancom_ASV-%j.err
#SBATCH -A p31629
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hannahschmitz2026@u.northwestern.edu

module purge all
module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0
module load gsl/2.7.1-gcc-12.3.0
module load fftw/3.3.10-gcc-12.3.0
module load gdal/3.7.0-gcc-12.3.0
module load nlopt/2.7.1-gcc-12.3.0
module load cmake/3.26.3-gcc-12.3.0

Rscript /gpfs/projects/b1052/hannah/r/scripts/run_ancombc2_ASV.R
