# PyCisTopic
**(if you want to test it, please read entirely before starting)**

pycisTopic is a Python framework for analyzing single-cell chromatin accessibility data (like scATAC-seq). It uses cis-topic modeling (topic modeling applied to genomic regions) to identify groups of co-accessible regions—called cis-topics—that represent shared regulatory programs across cells.

Documentation can be checked at: 
https://pycistopic.readthedocs.io/en/latest/index.html

## Notes for installation (SCENIC+)
I recommend installing pyCisTopic via SCENIC+ to ensure compatibility.
SCENIC+ was installed from GitHub (https://github.com/aertslab/scenicplus), which created a conda environment with pyCisTopic, scanpy/anndata, and all the dependencies. 

Then launched: 
````
conda activate scenicplus
jupyter notebook
````

## Dataset and Pseudobulk Strategy

I used a 10x Multiome (RNA + ATAC) Seurat object.

File: `integrated_prdm1_180924.rds`

Contains RNA-based cell type annotations (predicted.id) & ATAC data.

Although pycisTopic operates on ATAC data, the idea is to use RNA-derived annotations only to group cells in early steps.

## Why pseudobulk is required

scATAC-seq data is very sparse. To perform reliable peak calling, pycisTopic requires grouping cells and aggregating reads into pseudobulk profiles.

The thing is that peak calling cannot be done per single cell, we must first group cells → aggregate signal → call peaks

## Grouping strategies

pycisTopic does not require prior biological annotation, but it does require grouping. This can be done in two ways:

### 1. Supervised (the one we are trying to do)
Use existing annotations (RNA-based predicted.id)

### 2. Unsupervised
Cluster cells based only on chromatin accessibility - Use these clusters to define pseudobulks - Peaks are then called per cluster

## Conceptual workflow
1) Group cells (annotation or clustering)
2) Generate pseudobulk profiles per group
3) Perform consensus peak calling
4) Run pycisTopic at single-cell level using these peaks

**Important to undestand that:**

Pseudobulk is a technical requirement, not a final biological interpretation.
Downstream analysis (topics, motifs) remains single-cell.

## Cell Metadata Preparation (preprocess_pycistopic.R)

Before running pycisTopic, we generate a `cell_data.tsv` file that links cell barcodes to cell types and sample identity.

What this script does:
- Loads a Seurat object (multiome / integrated dataset)
- Subsets cells of interest (control cells with prefix ctl_)
- Removes Seurat-specific prefixes from barcodes
- Matches barcode format to fragment file (AAAC...-1)
- Creates a cell_data.tsv file required by pycisTopic

## Running pycisTopic (try_pycistopic.ipynb)

This notebook contains an end-to-end example of running pycisTopic, including both:

- Pseudobulk generation and peak calling
- Topic modeling on single-cell data

1. Pseudobulk generation (from fragments)

Inputs:

- `cell_data.tsv` (generated in preprocessing step)
- `atac_fragments.tsv.gz`
- `mm9.chrom.sizes`

Main function:
````
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
````
This step:

- Group cells according to cell_data.tsv
- Aggregates reads into pseudobulk profiles
- Prepares input for consensus peak calling

⚠️ Important:

This step is memory-intensive, and it may require running on a cluster (SLURM / HPC). I failed to run it locally. 
