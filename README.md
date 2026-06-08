16S rRNA Sequencing of Aerobic Granular Sludge, grouped by size.

Raw paired-end reads are available through the NCBI SRA, accession PRJNA1367764.

### 1.) Generate phyloseq object

Use define_ps.R 

### 2.) Determine significant taxa, using differential abundance analysis

Use run_ancombc2.R 

### 3.) Core Functions

- 01_load_ps.R: loads phyloseq object and defines necessary variables
- 02_metab_and_DA.R: defines functions for loading metabolism and differential abundance
  
### 4.) Plotting and Statistics
