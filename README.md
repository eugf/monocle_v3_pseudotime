# MONOCLE V3

### **REF**: [https://cole-trapnell-lab.github.io/monocle3/](https://cole-trapnell-lab.github.io/monocle3/)

### **PUB**: (Monocle v1) [https://www.nature.com/articles/nbt.2859](https://www.nature.com/articles/nbt.2859)

### **GOAL**: Run pseudotime analysis on single-cell RNA sequencing (scRNA-seq) data

---

FILE: monocle3-argparse-sbatch.sh

VARIABLES:


| VARIABLE        | EXAMPLE                                         | DESCRIPTION                                                                                                                                                               |
| :---------------- | :------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| SCRIPT_PATH     | ='../monocle_v3_pseudotime/monocle3-argparse.R' | File path to the`monocle3-argparse.R` script                                                                                                                              |
| VARIABLE        | ='age'                                          | The covariate we're intested in testing. Pseudotime will use the minimum of this column to to establish the root node to map the trajectories                             |
| DATA_TYPE       | ='rna'                                          | Data type. Choose 1 of the following: 'atac', 'rna'                                                                                                                       |
| CELL_TYPE       | ='VC'                                           | Cell type. Choose 1 of the following: 'Astro', 'ExN', 'InN', 'MG', 'Oligo', 'OPC', 'VC'                                                                                   |
| INPUT_PATH      | ='/INPUTS'                                      | Folder path where you will provide inputs for Monocle. This includes: a copy of the`.h5ad` file and files extracted from it to create Monocle's cell dataset (CDS) object |
| OUTPUT_PATH     | ='/OUTPUTS'                                     | Folder path where you intend to save Monocle's outputs                                                                                                                    |
| ANNDATA_PATH    | ='../01_anndata_object.h5ad'                    | File path to your anndata object saved in`.h5ad` format                                                                                                                   |
| GENES           | ='CLDN5,COLEC12,EPAS1,VCAM1'                    | String of gene markers separated by commas. If a gene is not found, it won't be plotted and no error will occur.                                                          |
| ALIGNMENT_GROUP | ='seq'                                          | Batch correction. Column must be part of your anndata file                                                                                                                |

# #* TODO - Make a README.md for each project with things that you might need to change and save which aren't in your scripts

explain how to use the Monocle script in my README.md once they've cloned
explain the variables
