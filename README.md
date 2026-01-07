# Lake_biodiversity_sequencingData_processing

This repository contains the bioinformatics pipeline and data processing scripts for environmental DNA (eDNA) metabarcoding analysis of freshwater eukaryotic communities. The raw sequencing data and metadata are deposited in NCBI under BioProject: 

# Project Overview
This pipeline processes Illumina MiniSeq paired-end sequencing data (150 bp) generated from water samples collected at the CNRS-ENS artificial lake platform in Nemours, France. The workflow amplifies two 18S rDNA regions (V7 and V9) using primer pairs Euka01 and Euka03 to assess multitrophic eukaryotic biodiversity.

#Repository Structure
.
├── OBITool4_scripts/          # OBITools 4 processing scripts
├── ngsfilter/                 # NGS filter files for demultiplexing
├── Obitools_output/           # Output files from OBITools 4 processing
├── MetabaR_Rcodes/            # R scripts for data cleaning with metabaR
├── PCR_GPS_files/             # PCR metadata files for metabaR
├── ngsfilter_files/           # Additional demultiplexing files for metabaR
├── sample_files/              # Sample metadata files for metabaR
├── Metabar_outputfiles/       # Cleaned output files from metabaR processing
└── 3-Merge_Euka0103_Replicates_DNA.R  # Script to merge Euka01 and Euka03 datasets

# Pipeline Workflow
## 1. Sequence Processing with OBITools 4
Raw sequencing data were processed using OBITools 4 following standard metabarcoding workflows.

### Scripts and files:

OBITool4_scripts/: Contains all OBITools 4 processing scripts
ngsfilter/: NGS filter files for sample demultiplexing
Obitools_output/: Processed sequence outputs

## 2. Data Cleaning with metabaR
Sequence data were cleaned and filtered using the metabaR R package to remove contaminants and low-quality sequences.

### Required files:

MetabaR_Rcodes/: R scripts for data cleaning pipeline
PCR_GPS_files/: PCR metadata
ngsfilter_files/: Demultiplexing information
sample_files/: Sample metadata

### Output:

Metabar_outputfiles/: Cleaned sequence data

## 3. Dataset Integration and Taxonomic Verification
Merging primer datasets:

### Script: 3-Merge_Euka0103_Replicates_DNA.R
Combines Euka01 (V7 region) and Euka03 (V9 region) datasets

### Taxonomic verification:
Each detected taxon was manually verified using:

GBIF (Global Biodiversity Information Facility)
WoRMS (World Register of Marine Species)
Published literature

Taxa were classified as aquatic, terrestrial, or both, and assessed for ecological plausibility in the study region.

# Final output:

1.Euka_DNA_Taxo_verified_pres_abs.csv: Taxonomically verified presence-absence matrix

# Requirements
Software

OBITools 4
R (≥ 4.3.3)

R Packages

metabaR (v1.0.0)
Additional dependencies as specified in individual R scripts

# Usage

Clone this repository
Process raw sequences using scripts in OBITool4_scripts/
Clean data using scripts in MetabaR_Rcodes/
Merge datasets using 3-Merge_Euka0103_Replicates_DNA.R
Final verified dataset: 1.Euka_DNA_Taxo_verified_pres_abs.csv

# Citation
If you use this pipeline or data, please cite:
[Your paper citation will go here once published]
Raw sequencing data: NCBI BioProject PRJNAXXXXXX

# License
[Choose an appropriate license - e.g., MIT, GPL-3.0, CC-BY-4.0]

# Contact
julie.guenat@unil.ch

# Acknowledgments
This work was conducted at the CNRS-ENS artificial lake platform in Nemours, France.
