source("./code/utils.R")
source("./code/downstream.R")
source("./code/normalize.R")
source("./code/hrmf_em.R")
source("./code/hrmf_init.R")
source("./code/image_tools.R")
#' Add edge score for 
#
# This function calculates an edge score for each spot based on emission probabilities
# and neighborhood information, classifies spots into "edge", "core", or normal, 
# and adds the results to the Seurat object.
#
# @param sample_dir Directory of the sample data (contains `seurat_object.rds`).
#   - `seurat_object.rds`: A Seurat object.
# @param spacna_dir Directory containing the output of spacna's result files (`cns.rds`, `count_norm.rds`, `bk_bic.rds`)
# @param plot_dir Directory to save plots and results
# @param tumor_content_dir Directory containing the output of tumor_content's result files (`tumor_content_df.rds`)
# @param neigh_cell List of neighbors for each spot (output of `get_spot_neighbors`).
# @param beta Numeric. Weight for Potts model in edge score calculation. Default is 1.
# @param edge_thre Numeric. Threshold for edge classification. Default is 0.25.
# @param tumor_thre Numeric. Threshold for tumor fraction (`tumor_content`). Default is -1.
# @return Seurat object with new metadata columns: `edge_score` and `edge`.
# @export
estimate_tumor_edge<- function(sample_dir, spacna_dir, plot_dir, tumor_content_dir,
                           beta = 1, edge_thre = 0.25, tumor_thre = -1) {
  obj<-readRDS(paste0(sample_dir,"seurat_object.rds"))
  cns <- readRDS(paste0(spacna_dir, "cns.rds"))
  count_norm <- readRDS(paste0(spacna_dir, "count_norm.rds"))
  gaussian_result<-readRDS(paste0(spacna_dir,'hrmf_para.rds'))
  df<-readRDS(paste0(tumor_content_dir,'tumor_content_df.rds'))

  obj <- AddMetaData(obj, df$p_estimate, col.name="tumor_content")
  spot_location <- get_spot_location(obj)
  neigh_cell <- get_spot_neighbors(spot_location)
  mus <- gaussian_result[[1]]
  sigmas <- gaussian_result[[2]]
  
  # Calculate emission probability matrix
  emission_prob_mat <- array(NA, dim = c(nrow(count_norm), ncol(count_norm), 3))
  for(j in 1:ncol(count_norm)){
    for(k in 1:3){
      emission_prob_mat[, j, k] <- dnorm(count_norm[,j], mus, sigmas)
    }
  }
  
  # Calculate Potts model probability matrix
  potts_prob_mat <- array(NA, dim = c(nrow(count_norm), ncol(count_norm), 3))
  for(i in 1:nrow(count_norm)){
    for(j in 1:ncol(count_norm)){
      tmp <- cns[i, neigh_cell[[j]]]
      for(k in 1:3){
        potts_prob_mat[i, j, k] <- sum(abs(tmp - k/2) < 0.05)
      }
    }
  }
  
  # Combine emission and Potts probabilities
  prob <- emission_prob_mat * exp(beta * potts_prob_mat)
  
  # Calculate entropy
  func <- function(x){
    x <- x / sum(x)
    return(-sum(x * log(x + 1e-6)))
  }
  entropy <- apply(prob, c(1,2), func)
  edge_score <- colMeans(entropy, na.rm = TRUE)
  
  # Add edge_score to Seurat object
  obj <- AddMetaData(obj, edge_score, col.name = "edge_score")
  p <- SpatialFeaturePlot(obj, features = "edge_score")
  ggplot2::ggsave(filename = file.path(plot_dir, "spacna_edge_score.png"), plot = p, width = 6, height = 5)
  
  # Classify spots
  edge <- rep("normal", ncol(obj))
  edge[obj$tumor_content > tumor_thre & obj$edge_score > edge_thre] <- "edge"
  edge[obj$tumor_content > tumor_thre & obj$edge_score <= edge_thre] <- "core"
  
  obj <- AddMetaData(obj, edge, col.name = "edge")
  p <- SpatialDimPlot(obj, "edge")
  ggplot2::ggsave(filename = file.path(plot_dir, "spacna_edge.png"), plot = p, width = 6, height = 5)

  
  return(obj)
}

# sample_dir <-  ""
# plot_dir<- ""
# spacna_dir <- ""

# tumor_content<-estimate_tumor_edge(sample_dir =sample_dir,plot_dir = plot_dir,spacna_dir = spacna_dir,tumor_content_dir=spacna_dir)