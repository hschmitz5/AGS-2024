16S rRNA Sequencing of Aerobic Granular Sludge, grouped by size.

Raw paired-end reads are available through the NCBI SRA, accession PRJNA1367764.

### 1.) Generate phyloseq object

- define_ps.R -- generates ASV and genus level data
- define_ps_metab -- generates functional level data

### 2.) Determine significant taxa, using differential abundance (DA) analysis

- run_ancombc2_genus.R -- runs DA at the genus level
- run_ancombc2_metab.R -- runs DA at the functional group level

### 3.) Core Functions

- 01_load_ps.R -- loads phyloseq object and defines necessary variables
- 02_sum_rel_ab_by_function.R -- agglomerates relative abundance by functional group
  
### 4.) Plotting and Statistics
