# MONOCLE V3

### REF: [https://cole-trapnell-lab.github.io/monocle3/](https://cole-trapnell-lab.github.io/monocle3/)

### PUB: (Monocle v1) [https://www.nature.com/articles/nbt.2859](https://www.nature.com/articles/nbt.2859)

### GOAL: Run pseudotime analysis on single-cell RNA sequencing (scRNA-seq) data

---

**FILE**: monocle3-argparse-sbatch.sh

**SBATCH PARAMETERS:**


| COMMAND           | DESCRIPTION                                                                                                                                       |
| ------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| --cpus-per-task=1 | Parts of Monocle are limited to 1 CPU. You can run faster at the expense of reproducibility because it will randomly subsample if you use >1 CPU. |
| --mem=600g        | 500,000 cells requires ~600 GB of RAM to run                                                                                                      |
| --time=24:00:00   | Time and memory scales with the number of principal components used as well                                                                       |

**VARIABLES:**


| VARIABLE        | EXAMPLE                                         | DESCRIPTION                                                                                                                                                               |
| :---------------- | :------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| SCRIPT_PATH     | ='../monocle_v3_pseudotime/monocle3-argparse.R' | File path to the`monocle3-argparse.R` script                                                                                                                              |
| VARIABLE        | ='age'                                          | The covariate we're interested in testing. Pseudotime will use the minimum of this column to to establish the root node to map the trajectories                           |
| DATA_TYPE       | ='rna'                                          | Data type. Choose 1 of the following: 'atac', 'rna'                                                                                                                       |
| CELL_TYPE       | ='VC'                                           | Cell type. Choose 1 of the following: 'Astro', 'ExN', 'InN', 'MG', 'Oligo', 'OPC', 'VC'                                                                                   |
| INPUT_PATH      | ='/INPUTS'                                      | Folder path where you will provide inputs for Monocle. This includes: a copy of the`.h5ad` file and files extracted from it to create Monocle's cell dataset (CDS) object |
| OUTPUT_PATH     | ='/OUTPUTS'                                     | Folder path where you intend to save Monocle's outputs                                                                                                                    |
| ANNDATA_PATH    | ='../01_anndata_object.h5ad'                    | File path to your anndata object saved in`.h5ad` format                                                                                                                   |
| GENES           | ='CLDN5,COLEC12,EPAS1,VCAM1'                    | String of gene markers separated by commas. If a gene is not found, it won't be plotted and no error will occur.                                                          |
| ALIGNMENT_GROUP | ='seq'                                          | Batch correction. Column must be part of your anndata file                                                                                                                |

### NOTES:

- On RStudio, Monocle 3 was limited to R/4.1
- As of 8/29/24, Monocle 3 will run on up to R/4.4 when running thru the command line

# HOW TO RUN:

1) Fill in the sbatch parameters that you think are reasonable for the size of your dataset
2) Fill in the variables in the `monocle2-argparse-sbatch.sh` script file with the paths and variables you want to test with your data
3) Open a terminal and navigate to where you've stored the `monocle2-argparse-sbatch.sh` script file
4) Run using:

```
sbatch monocle2-argparse-sbatch.sh
```

5) The `/OUTPUTS` folder will have: