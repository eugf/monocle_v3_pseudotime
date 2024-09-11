#!/usr/bin/env Rscript

# IMPORTS
require(argparse)
require(ggplot2)
require(monocle3)
require(rhdf5)
require(remotes)
if (!require(scrattch.io)) {
    remotes::install_github("AllenInstitute/scrattch.io", upgrade=TRUE)
}
library(scrattch.io)
require(data.table)

set.seed(1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# INITIALIZE ARGPARSE
parser <- ArgumentParser()

# DEFINE ARGUMENTS THE USER WILL PROVIDE
parser$add_argument(
    '--alignment-group', 
    type        = 'character',
    default     = 'seq'
)

parser$add_argument(
    '--variable', 
    help        = 'NOTE: this is not really being used right now, but it is for the experimental variable we want to look at (ex: PMI, Age) and will sort the pseudotime by',
    type        = 'character', 
    default     = 'age'
)

parser$add_argument(
    '--data-type', 
    required    = TRUE,
    type        = 'character', 
    default     = 'rna',
    choices     = c('atac', 'rna')
)

parser$add_argument(
    '--cell-type', 
    required    = TRUE,
    type        = 'character', 
    choices     = c('Astro', 'ExN', 'InN', 'MG', 'Oligo', 'OPC', 'VC')
)

parser$add_argument(
    '--input-path', 
    help        = 'Base folder where anndata file and converted files should be stored to be used as inputs for downstream operations',
    required    = TRUE,
    type        = 'character'
)

parser$add_argument(
    '--output-path', 
    help        = 'Base folder for processed files and plots to be stored',
    required    = TRUE,
    type        = 'character'
)

parser$add_argument(
    '--anndata-path', 
    help        = "WARNING: currently using my copies of Shahroze's anndata files, but this is intended to take a full anndata file path and make a copy of it for my local input path to operate on",
    required    = TRUE,
    type        = 'character'
)

parser$add_argument(
    '--genes', 
    required    = TRUE,
    type        = 'character',
    nargs       = '+'
)

# LOAD USER-PROVIDED ARGUMENTS
args <- parser$parse_args()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# REMEMBER: the order matters on some of these!

# DEFINE VARIABLES
alignment_group <- args$alignment_group
num_dim       <- 50             # Default = 50, will filter this later
current_date  <- as.character(Sys.Date())
current_time  <- format(Sys.time(),"%H-%M-%S")

# INPUT PATH
data_type     <- args$data_type
cell_type     <- args$cell_type 
input_path    <- file.path(args$input_path, data_type, cell_type)

# OUTPUT PATH
column_of_interest  <- args$variable
output_path         <- file.path(args$output_path, column_of_interest, data_type, cell_type, current_date)

# ANNDATA PATH
h5_path             <- args$anndata_path
anndata_file_name   <- basename(h5_path)
anndata_split       <- strsplit(anndata_file_name, "\\.")[[1]]
anndata_name_only   <- anndata_split[1]
anndata_file_ext    <- anndata_split[2]

#* TODO - eventually may not want to make extra copies of the anndata files
# ANNDATA COPY
# CHECK - is the anndata in my own intended local directory for inputs already or somewhere else like someone's Biowulf folder? Make a copy if elsewhere
comparison_path     <- file.path(input_path, anndata_file_name)
if (h5_path == comparison_path) {
  cat('h5ad paths match! do nothing \n')
} else {
  dir.create(dirname(comparison_path), recursive = TRUE)
  file.copy(h5_path, comparison_path, overwrite = TRUE)
  cat('h5ad file copied to: ', comparison_path, '\n')
  h5_path <- comparison_path
}

# OUTPUT - CDS PATH
cds_path      <- file.path(output_path, paste(sep='-', "cds", anndata_name_only, cell_type, ".rds"))
cds_path_wip  <- file.path(output_path, paste(sep='-', "cds", anndata_name_only, cell_type, "wip.rds"))
cat('cds_path     = \n', cds_path, '\n')
cat('cds_path_wip = \n', cds_path_wip, '\n')

# MARKER GENE LIST PROCESSING
genes <- args$genes
cat('genes = ', unlist(list(genes)), '\n')

# Initialize an empty list to store genes
genes_list <- character()

# NOTE: R doesn't like the list provided from the bash script, so this will take that single string and turn it into a list in R
for (gene in strsplit(genes, ",")) {
  genes_list <- c(genes_list, gene)
}
cat(unlist(genes_list))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DEFINE FUNCTIONS 

# FILE NAMING - Construct the file name
my_plot_name <- function(plot_name = "") {

  # Construct file name and output folder based on user-input and other defined vars
  file_name = paste0('PLOT-', cell_type, '-', plot_name, '.png')
  
  # Construct the final file path
  file_path = file.path(output_path, file_name)
  cat('file_path = \n', file_path, '\n')

  # CHECK - if the output directory exists, if not: create it
  if (!dir.exists(dirname(file_path))) {
    dir.create(dirname(file_path), recursive = TRUE)
  }

  return(file_path)
}

# SAVE PLOTS
# - Include some sensible defaults for width, height, DPI
# - Take in user-defined name, save last plot
my_plot_save <- function(
    filename  = my_plot_name(),
    width     = 8, 
    height    = 6, 
    dpi       = 300
    ) 
  {  
  # Define the file format and dimensions to open with graphics device
  ggsave(
    filename  = filename,
    width     = width, 
    height    = height, 
    dpi       = dpi
    )
  }

# CONVERT - nested H5AD columns from the OBS slot into a flattened DF
flatten_df <- function(df) {

    column_names  <- names(df)
    indices       <- df$`_index`
    N             <- length(indices)

    # Initialize new DF based on the full length of input table
    out_df <- data.frame(barcode = indices)

    # Iterate thru each column name
    for(i in column_names) {

        mycols <- names(df[[i]])

        # If column needs to be flattened
        if(length(mycols) == 2) {
            vals <- df[[i]]$categories[(1 + df[[i]]$codes)]

        # If column is already flattened
        } else if(length(df[[i]]) == N) {
            vals <- df[[i]]

        # Just in case error
        } else {
            cat('ERROR: Vector inside OBS object is unexpected length \n')
            exit()
        }

        # Put flattened column into the DF
        out_df[[i]] <- vals
    }

    # Output final reformatted DF
    return(out_df)
}

# CONVERT - ANNDATA TO SPARSE MATRIX 
# Remake the h5ad read function from: 
# https://rdrr.io/github/AllenInstitute/scrattch.io/src/R/read_h5ad.R
read_h5ad_dgCMatrix <- function(h5ad, target='X') {

  #library(Matrix)

  # Open H5AD file
  root    <- rhdf5::H5Fopen(h5ad)

  # Extract metadata and data
  i_path  <- paste0(target,"/indices")
  p_path  <- paste0(target,"/indptr")
  x_path  <- paste0(target,"/data")

  cat("Reading indices \n")
  i <- read_tome_vector(root, i_path)
  cat("Reading pointers \n")
  p <- read_tome_vector(root, p_path)
  cat("Reading values \n")
  x <- read_tome_vector(root, x_path)

  cat("Reading observations \n")
  o <- rhdf5::h5read(root, "/obs")
  cat("Reading variables \n")
  v <- as.vector(rhdf5::h5read(root, "/var")$`_index`)

  cat("Reading dimensions \n")
  # dims <- c(length(v), length(o))

  #! WARNING - should close file or will have access locked issues
  H5Fclose(root)

  # FUNCTION - Flatten the OBS DF
  o2 <- flatten_df(o)

  rownames(o2) <- o$`_index`
  dims <- c(length(v), nrow(o2))

  cat("Assembling dgCMatrix \n")
  m <- Matrix::sparseMatrix(
    i = i,
    p = p,
    x = x,
    index1 = FALSE,
    dims = dims)
  rownames(m) <- v
  colnames(m) <- rownames(o2)

  # OUTPUTS: X, obs, and vars
  output <- list(v, o2, m)
  return(output)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RUN H5AD CONVERSION
# - Read in H5AD file as a sparse matrix
dat <- read_h5ad_dgCMatrix(h5_path, target='X')

# Need to access the cell_metadata file later
cell_metadata = dat[[2]]

#! WARNING - should name it 'gene_short_name'
gene_data <- data.frame(gene_short_name = dat[[1]])
rownames(gene_data) <- gene_data$gene_short_name

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CREATE MONOCLE CDS - using the outputs from the H5AD conversion
cds <- new_cell_data_set(
  as(dat[[3]], "sparseMatrix"),
  cell_metadata = cell_metadata,
  gene_metadata = gene_data)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CHECK - if the parent directory exists, if not: create it
if (!dir.exists(dirname(cds_path))) {
  dir.create(dirname(cds_path), recursive = TRUE)
}

# EXPORT - CDS file
saveRDS(cds, file = cds_path)
  
# CHECK - CDS.rds file path to save to
cat('Saved CDS to: \n', cds_path, '\n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MONOCLE WORKFLOW
cat('Start Monocle workflow \n')

# # RELOAD - CDS (if needed)
# cds <- readRDS(file = cds_path)

# PREPROCESS THE DATA
cds <- preprocess_cds(cds, num_dim = num_dim)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT - check that you're using enough PCs
plot_pc_variance_explained <- plot_pc_variance_explained(cds)
plot_pc_variance_explained +
  labs(title = paste(sep='-', cell_type, '-pc_variance_explained'))

# SAVE PLOT
save_plot_pc_variance_explained <- my_plot_name('pc_variance_explained')
my_plot_save(save_plot_pc_variance_explained)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CALCULATE - the top 95% PC contributors
# Extract just the final PCA values (it's already sorted in order of largest to smallest)
pca_list <- cds@reduce_dim_aux@listData[["PCA"]]@listData[["model"]]@listData[["prop_var_expl"]]

# Use cumulative sum to add up each PC's value
# Turn into a boolean if the 95% threshold has been met
total_var_ex <- cumsum(pca_list) >= 0.95

# Get the index of the 1st element that satisfies that condition
pca_list_95 <- which(total_var_ex)[1]

# CHECK
cat('pca_list = \n', pca_list, '\n')
cat('pca_list[pca_list_95] = \n', pca_list[pca_list_95], '\n')

# REDEFINE - set the number of dimensions to use that index for 95%
num_dim <- pca_list_95
cat('NEW num_dim = \n', num_dim, '\n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ALIGNMENT - correct for batch effects using the provided alignment group (default = sequencing group)
cds <- align_cds(cds, num_dim = num_dim, alignment_group = alignment_group)

# # REDUCE DIMENSIONALITY
# cds <- reduce_dimension(cds, reduction_method='UMAP')
#! IMPORTANT! - Use these settings to prevent randomness
cds <- reduce_dimension(cds, umap.fast_sgd = FALSE, cores = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#* TODO - loop this???
# # ENABLE - only for the big lists of gene markers
# counter <- 1

# for (item in genes_big_list) {
#   genes = item
#   # cat(genes, '\n')

#   # PLOT - by gene markers of interest (user-supplied list)
#   plot_gene_markers <- plot_cells(cds, genes = genes, show_trajectory_graph=FALSE) +
#     labs(title = paste(sep='-', cell_type, 'gene_markers', counter))

#   # SAVE PLOT
#   save_plot_gene_markers <- my_plot_name(paste0('gene_markers-', counter))
#   my_plot_save(save_plot_gene_markers)

#   # Increment the count by 1
#   counter <- counter + 1
#   # cat(counter, '\n')
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # PLOT - by gene markers of interest (user-provided list)
# plot_gene_markers <- plot_cells(cds, genes = genes)
plot_gene_markers <- plot_cells(cds, genes = genes_list, show_trajectory_graph = FALSE) +
  labs(title = paste(sep='-', cell_type, 'gene_markers'))

# SAVE PLOT
save_plot_gene_markers <- my_plot_name('gene_markers')
my_plot_save(save_plot_gene_markers)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CLUSTER CELLS (required! even if you're not subclustering!)
cds <- cluster_cells(cds, reduction_method='UMAP', k=20)

# #! Warning message:
# In cluster_cells_make_graph(data = data, weight = weight, cell_names = cell_names,  :
#   The nearest neighbors includes the point itself, k must be smaller than
# the total number of points - 1 (all other points) - 1 (itself)! Total number of points is 26

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT - by Age
plot_age <- plot_cells(cds, color_cells_by = column_of_interest) +
  labs(title = paste(sep='-', cell_type, column_of_interest))

# SAVE PLOT
save_plot_age <- my_plot_name(column_of_interest)
my_plot_save(save_plot_age)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # TEST
# # PLOT - by column of interest (that's the covariates we are interested in looking at and that pseudotime will filter by)
# plot_age <- plot_cells(cds, color_cells_by = column_of_interest) +
#   labs(title = paste(sep='-', cell_type, column_of_interest))

# # SAVE PLOT
# save_plot_age <- my_plot_name(column_of_interest)
# my_plot_save(save_plot_age)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT - by Brain Bank
plot_brain_bank <- plot_cells(
  cds, 
  color_cells_by = "brain_bank",
  label_cell_groups = FALSE) +
  labs(title = paste(sep='-', cell_type, 'brain_bank'))

# SAVE PLOT
save_plot_brain_bank <- my_plot_name('brain_bank')
my_plot_save(save_plot_brain_bank)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # TEST - enable for debugging if needed
# # EXPORT CDS (save point before starting pseudotime)
# saveRDS(cds, file = cds_path_wip)
# # RELOAD - CDS (if needed)
# cds <- readRDS(file = cds_path_wip)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PSEUDOTIME
cat('Start pseudotime workflow \n')

# Learn the principal graph
cds <- learn_graph(cds)

# # Order cells in pseudotime based on the covariates column you're interested in
# # random_root <- sample(colnames(cds), 1)
# roots <- rownames(cell_metadata[cell_metadata$age == min(cell_metadata$age), ])
roots <- rownames(cell_metadata[cell_metadata[[column_of_interest]] == min(cell_metadata[[column_of_interest]]), ])

cds <- order_cells(cds, root_cells=roots)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT - pseudotime trajectory, but reduce the #s
plot_pseudotime_reduce_numbers <- plot_cells(
  cds,
  color_cells_by          = "pseudotime",
  label_cell_groups       = FALSE,
  label_leaves            = FALSE,
  label_branch_points     = FALSE,
  graph_label_size        = 1.5,
  label_roots             = TRUE,
  label_groups_by_cluster = TRUE
  # label_principal_points  = TRUE
  # show_trajectory         = FALSE
) +
  labs(title = paste(sep='-', cell_type, 'pseudotime_reduce_numbers'))

# SAVE PLOT
save_plot_pseudotime_reduce_numbers  <- my_plot_name('pseudotime_reduce_numbers')
my_plot_save(save_plot_pseudotime_reduce_numbers)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT - pseudotime trajectory, but remove the #s entirely
plot_pseudotime_wo_nos <- plot_cells(
  cds,
  color_cells_by      = "pseudotime",
  label_cell_groups   = FALSE,
  label_leaves        = FALSE,
  label_branch_points = FALSE,
  graph_label_size    = 1.5,
  show_trajectory     = FALSE
) +
  labs(title = paste(sep='-', cell_type, 'pseudotime_wo_nos'))

# SAVE PLOT
save_plot_pseudotime_wo_nos <- my_plot_name('pseudotime_wo_nos')
my_plot_save(save_plot_pseudotime_wo_nos)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT - pseudotime trajectory
plot_pseudotime <- plot_cells(cds, color_cells_by = 'pseudotime') +
  labs(title = paste(sep='-', cell_type, 'pseudotime'))

# SAVE PLOT
save_plot_pseudotime <- my_plot_name('pseudotime')
my_plot_save(save_plot_pseudotime)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NOTE - Shahroze used this for plotting, saves it as a col in the metadata too
cell_metadata$pseudotime <- pseudotime(cds)

# EXPORT - updated cell_metadata as a TSV file
fwrite(
  cell_metadata,
  file      = file.path(output_path, 'cell_metadata.tsv'), 
  sep       = '\t', 
  row.names = FALSE
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # PLOT
# plot_pseudotime_histogram_freq_age <- 
#   ggplot(cell_metadata, aes(x=age, y=pseudotime)) + 
#   geom_point() + 
#   theme_minimal() + 
#   labs(title = 'Pseudotime VS Age',
#        x     = 'Age',
#        y     = 'Pseudotime') +
#   theme(panel.background = element_rect(fill = "white"))

# # SAVE PLOT
# save_plot_pseudotime_histogram_freq_age <- my_plot_name('pseudotime_histogram_freq_age')
# my_plot_save(save_plot_pseudotime_histogram_freq_age)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
plot_pseudotime_histogram_freq_age_jitter <- 
  ggplot(cell_metadata, aes(x=age, y=pseudotime)) + 
  geom_jitter(alpha=0.7) + 
  theme_minimal() + 
  labs(title = 'Pseudotime VS Age',
       x     = 'Age',
       y     = 'Pseudotime') +
  theme(panel.background = element_rect(fill = "white"))

# SAVE PLOT
save_plot_pseudotime_histogram_freq_age_jitter <- my_plot_name('pseudotime_histogram_freq_age_jitter')
my_plot_save(save_plot_pseudotime_histogram_freq_age_jitter)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
plot_pseudotime_histogram_dist <- 
  ggplot(cell_metadata, aes(x=pseudotime)) + 
  geom_histogram(binwidth=0.5, fill='steelblue', color='black', alpha= 0.7) +
  theme_minimal() + 
  labs(title = 'Distribution of Pseudotime',
       x       = 'Pseudotime',
       y       = 'Frequency')

# SAVE PLOT
save_plot_pseudotime_histogram_dist <- my_plot_name('pseudotime_histogram_dist')
my_plot_save(save_plot_pseudotime_histogram_dist)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOGGING 
# - save job metadata, user-provided inputs, and important calculated variables to a plaintext log file

# CREATE - filename + path where you want to save the log
my_log_file <- file.path(output_path, paste0('log-', current_date, '_', current_time, '.txt'))

# OPEN SINK (for logging terminal outputs to a text file)
sink(my_log_file)

# JOB INFO:
script_name <- deparse(substitute(sys.call(-1)))
cat("JOB INFO: \n")
cat(paste0("current_datetime  = ", current_date, "_", current_time, " \n"))
cat("script_name       = ", script_name, " \n")

# ARGPARSE USER INPUTS
cat('\nUSER-PROVIDED ARGUMENTS: \n')
print(args)

# VARIABLES CALCULATED IN SCRIPT
cat('\nVARIABLES CALCULATED IN SCRIPT: \n')
cat('95% of PCs  = \n', num_dim, ' PCs used', '\n')
cat('Pseudotime ordered by minimum', column_of_interest, '= \n', min(cell_metadata[[column_of_interest]]), '\n')

# FINAL OUTPUT FOLDER PATH
cat('\nOUTPUTS SAVED TO: \n', output_path, '\n')

# CLOSE SINK
sink(NULL)