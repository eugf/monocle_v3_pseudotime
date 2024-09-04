#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=600g 
#SBATCH --partition=largemem
#SBATCH --mail-type=BEGIN,END 
#SBATCH --time=24:00:00

#! WARNING - 500K cells requires ~600 GB of RAM

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DEFINE VARIABLES FOR ARGPARSE - change these!!!
SCRIPT_PATH='/data/CARD_singlecell/MONOCLE_V3/SRC/MVP/FINAL/monocle_v3_pseudotime/monocle3-argparse.R'
VARIABLE='age'
DATA_TYPE='rna'
CELL_TYPE='VC'
INPUT_PATH='/data/CARD_singlecell/MONOCLE_V3/INPUTS/PFC_brain_atlas'
OUTPUT_PATH='/data/CARD_singlecell/MONOCLE_V3/OUTPUTS/PFC'
ANNDATA_PATH='/data/CARD_singlecell/MONOCLE_V3/INPUTS/PFC_brain_atlas/rna/VC/01_anndata_object.h5ad'
GENES='CLDN5,COLEC12,EPAS1,VCAM1'
ALIGNMENT_GROUP='seq'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Output script name and the initial command-line arguments to slurm file 
echo "SCRIPT NAME   : $0"
echo "JOB INFO      : $@"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NOTES: 
# - previously on RStudio, Monocle 3 was limited to R/4.1
# - currently (8/29/24) Monocle 3 will run on up to R/4.4

# IMPORTS
module load R/4.4

# RUN THIS SCRIPT
Rscript \
    ${SCRIPT_PATH} \
    --variable ${VARIABLE} \
    --data-type ${DATA_TYPE} \
    --cell-type ${CELL_TYPE} \
    --input-path ${INPUT_PATH} \
    --output-path ${OUTPUT_PATH} \
    --anndata-path ${ANNDATA_PATH} \
    --genes ${GENES} \
    --alignment-group ${ALIGNMENT_GROUP} || (echo 'R failed!'; exit 1)