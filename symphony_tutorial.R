install.packages("symphony")

devtools::install_github("immunogenomics/symphony")

library(symphony)
suppressPackageStartupMessages({
  # Analysis
  library(harmony)
  library(irlba)
  library(data.table)
  library(dplyr)
  
  # Plotting
  library(ggplot2)
  library(ggthemes)
  library(ggrastr)
  library(RColorBrewer)
  })

plotBasic = function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                     title = 'Query',         # Plot title
                     color.by = 'cell_type',  # metadata column name for coloring
                     facet.by = NULL,         # (optional) metadata column name for faceting
                     color.mapping = NULL,    # custom color mapping
                     legend.position = 'right') {  # Show cell type legend
  
  p = umap_labels %>%
    dplyr::sample_frac(1L) %>% # permute rows randomly
    ggplot(aes(x = UMAP1, y = UMAP2)) + 
    geom_point_rast(aes(col = get(color.by)), size = 1, stroke = 0.4, shape = 16)
  if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
  
  # Default formatting
  p = p + theme_bw() +
    labs(title = title, color = color.by) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position=legend.position) +
    theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')
  
  if(!is.null(facet.by)) {
    p = p + facet_wrap(~get(facet.by)) +
      theme(strip.text.x = element_text(size = 12)) }    
  return(p)
}

# Colors for PBMCs
pbmc_colors = c("B" = "#66C2A5", "DC" = "#FC8D62", "HSC" = "#8DA0CB", "MK" = "#E78AC3", 
                "Mono_CD14" = "#A6D854", "Mono_CD16" = "#f2ec72", "NK" = "#62AAEA", 
                "T_CD4" = "#D1C656", "T_CD8" = "#968763")

#load the data
load('./data/pbmcs_exprs_small.rda')
load('./data/pbmcs_meta_small.rda')

dim(pbmcs_exprs_small)
dim(pbmcs_meta_small)

pbmcs_meta_small %>% head(4)

idx_query = which(pbmcs_meta_small$donor == "5p") # use 5' dataset as the query
ref_exp_full = pbmcs_exprs_small[, -idx_query]
ref_metadata = pbmcs_meta_small[-idx_query, ]
query_exp = pbmcs_exprs_small[, idx_query]
query_metadata = pbmcs_meta_small[idx_query, ]

# Build Symphony reference

## Option 1: Build from Harmony object (preferred method)

ref_exp_full[1:5, 1:2] # Sparse matrix with the normalized genes x cells matrix

var_genes = vargenes_vst(ref_exp_full, groups = as.character(ref_metadata[['donor']]), topn = 1000)
ref_exp = ref_exp_full[var_genes, ]
dim(ref_exp)

vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
vargenes_means_sds$stddev = rowSDs(ref_exp, vargenes_means_sds$mean)
head(vargenes_means_sds)

ref_exp_scaled = scaleDataWithStats(ref_exp, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)

set.seed(0)
s = irlba(ref_exp_scaled, nv = 20)
Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
loadings = s$u

set.seed(0)
ref_harmObj = harmony::HarmonyMatrix(
  data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
  meta_data = ref_metadata, ## dataframe with cell labels
  theta = c(2),             ## cluster diversity enforcement
  vars_use = c('donor'),    ## variable to integrate out
  nclust = 100,             ## number of clusters in Harmony model
  max.iter.harmony = 20,    ## max number of iterations
  return_object = TRUE,     ## return the full Harmony model object
  do_pca = FALSE            ## don't recompute PCs
)  

reference = buildReferenceFromHarmonyObj(
  ref_harmObj,            # output object from HarmonyMatrix()
  ref_metadata,           # reference cell metadata
  vargenes_means_sds,     # gene names, means, and std devs for scaling
  loadings,               # genes x PCs matrix
  verbose = TRUE,         # verbose output
  do_umap = TRUE,         # set to TRUE to run UMAP
  save_uwot_path = './testing_uwot_model_1') # file path to save uwot model