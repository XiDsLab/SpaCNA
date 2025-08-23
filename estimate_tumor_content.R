source("./code/utils.R")
source("./code/downstream.R")
source("./code/normalize.R")
source("./code/hrmf_em.R")
source("./code/hrmf_init.R")
source("./code/image_tools.R")
# Estimate tumor content using Cell2location and SpaCNA
# @param sample_dir Directory of the sample data (contains `seurat_object.rds`).
#   - `seurat_object.rds`: A Seurat object.
# @param spacna_dir Directory containing the output of spacna's result files (`cns.rds`, `count_norm.rds`, `bk_bic.rds`)
# @param plot_dir Directory to save plots and results
# @param K Number of clones for clustering (default: 7)
# @param cnv_thre Minimum number of CNA segments to consider a clone as tumor (default: 10)
# @return Seurat object with metadata `tumor_content` added
estimate_tumor_content <- function(sample_dir, spacna_dir, plot_dir, K = 7, cnv_thre = 10) {
  obj<-readRDS(paste0(sample_dir,"seurat_object.rds"))
  if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  cns <- readRDS(file.path(spacna_dir, "cns.rds"))
  count_norm <- readRDS(file.path(spacna_dir, "count_norm.rds"))
  bk_bic <- readRDS(file.path(spacna_dir, "bk_bic.rds"))
  cell_barcodes <- colnames(count_norm)
  
  temp <- count_norm
  d <- parallelDist::parallelDist(t(temp))
  hc_re <- hclust(d, method="ward.D2")
  clone <- cutree(hc_re, k=K)
  
  plot_heatmap(t(cns),
               cell_cluster=clone,
               cell_annotation=clone,
               output_dir=plot_dir,
               output_name="cns_clone_annot.png")
  
  n_seg <- length(bk_bic) - 1
  cns_seg <- matrix(nrow = ncol(cns), ncol = n_seg)
  rownames(cns_seg) <- colnames(cns)
  for(i in 1:n_seg){
    cns_seg[,i] <- cns[bk_bic[i],]
  }
  
  copy_ratio_seg <- cns_seg
  for(i in 1:ncol(copy_ratio_seg)){
    copy_ratio_seg[,i] <- apply(count_norm[bk_bic[i]:bk_bic[i+1],], 2, mean)
  }
  
  df <- data.frame(
    p_estimate = rep(0, ncol(cns)),
    p_estimate_confidence = rep(0, ncol(cns)),
    row.names = colnames(cns)
)  
  copy_ratio_pure_gain <- rep(1, n_seg)
  copy_ratio_pure_loss <- rep(1, n_seg)
  for(i in 1:n_seg){
    tmp <- copy_ratio_seg[, i]
    f <- cns_seg[,i] < 1
    if(sum(f) > 50){
      copy_ratio_pure_loss[i] <- mean(sort(tmp[f])[1:10])
    }
    f <- cns_seg[,i] > 1
    if(sum(f) > 50){
      copy_ratio_pure_gain[i] <- mean(sort(tmp, decreasing = T)[1:10])
    }
  }
  
  get_cnv_clone <- function(x, thre = 0.3){
    if(sum(x == 0.5)/length(x) > thre){
      return(0.5)
    } else if(sum(x == 1.5)/length(x) > thre){
      return(1.5)
    } else {
      return(1)
    }
  }
  
  for(k in seq_len(K)){
    clone_cnv <- apply(cns_seg[clone==k,], 2, get_cnv_clone)
    if(sum(clone_cnv != 1) < cnv_thre){
      next
    }
    copy_ratio_pure <- rep(1, length(clone_cnv))
    f <- clone_cnv > 1
    copy_ratio_pure[f] <- copy_ratio_pure_gain[f]
    f <- clone_cnv < 1
    copy_ratio_pure[f] <- copy_ratio_pure_loss[f]
    
    clone_cell <- cell_barcodes[clone == k]
    for(i in seq_along(clone_cell)){
      f <- clone_cnv != 1
      df_ <- data.frame(x = copy_ratio_pure[f] - 1, y = copy_ratio_seg[clone_cell[i], f] - 1)
      re <- lm(y ~ x + 0, data = df_)
      df[clone_cell[i], "p_estimate"] <- re$coefficients
      df[clone_cell[i], "p_estimate_confidence"] <- summary(re)$adj.r.squared
    }
  }
  
  ggplot(df, aes(x=p, y=p_estimate)) + geom_point()
  
  obj <- AddMetaData(obj, df$p_estimate, col.name="tumor_content")
  p <- SpatialFeaturePlot(obj, features = "tumor_content")
  ggplot2::ggsave(filename = file.path(plot_dir, "spacna_tumor_content.png"), plot = p, width = 6, height = 5)
  
  # Save df
  saveRDS(df, file.path(plot_dir, "tumor_content_df.rds"))
  return(obj)
}


# sample_dir <-  ""
# plot_dir<- ""
# spacna_dir <- ""

# tumor_content<-estimate_tumor_content(sample_dir =sample_dir,plot_dir = plot_dir,spacna_dir = spacna_dir,tumor_content_dir=spacna_dir)