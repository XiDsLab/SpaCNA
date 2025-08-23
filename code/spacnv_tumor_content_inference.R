source("./code_ver2/utils.R")
source("./code_ver2/downstream.R")
source("./code_ver2/normalize.R")
source("./code_ver2/hrmf_em.R")
source("./code_ver2/hrmf_init.R")

sample_dir <- "./data/CRC_2/"
plot_dir <- "./CRC_2_analysis/"
obj <- readRDS("./data/CRC_2/seurat_object.rds")

## input
cns <- readRDS(paste0(sample_dir, "cns.rds")) # copy number state
count_norm <- readRDS(paste0(sample_dir, "count_norm.rds")) # copy ratio
bk_bic <- readRDS(paste0(sample_dir, "bk_bic.rds")) # breakpoints
cell_barcodes <- colnames(count_norm)
K <- 7 # number of cluster
cnv_thre <- 10 # number of minimum cnv segment 

## cluster
temp <- count_norm
d <- parallelDist::parallelDist(t(temp))
hc_re <- hclust(d, method="ward.D2")
clone <- cutree(hc_re, k=K)

obj <- AddMetaData(obj, clone, col.name="clone")
SpatialDimPlot(obj, group="clone")

plot_heatmap(t(cns),
             cell_cluster=clone,
             cell_annotation=clone,
             output_dir=plot_dir,
             output_name="cns_clone_annot_ver2.png")

## merge cell-bin copy number state and copy ratio into cell-segment matrix 
n_seg <- length(bk_bic) - 1
cns_seg <- matrix(nrow = ncol(cns),
                  ncol = n_seg)
rownames(cns_seg) <- colnames(cns)
for(i in 1:n_seg){
  cns_seg[,i] <- cns[bk_bic[i],]
}

copy_ratio_seg <- cns_seg
for(i in 1:ncol(copy_ratio_seg)){
  copy_ratio_seg[,i] <- apply(count_norm[bk_bic[i]:bk_bic[i+1],], 2, mean)
}

## tumor content inference 
tumor_content <- data.frame(row.names = colnames(cns),
                            p <- rep(0, ncol(cns)))
df <- tumor_content
df$p_estimate <- 0
df$p_estimate_confidence <- 0 # adjusted R^2 in linear regression

# copy ratio baseline
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

# function to decide common CNV of a cell cluster
get_cnv_clone <- function(x, thre = 0.5){
  if(sum(x == 0.5)/length(x) > thre){
    return(0.5)
  }else if(sum(x == 1.5)/length(x) > thre){
    return(1.5)
  }else{
    return(1)
  }
}

# linear regression
for(k in c(1:K)){
  # print(sum(clone==k))
  clone_cnv <- apply(cns_seg[clone==k,], 2, get_cnv_clone)
  print(sum(clone_cnv!=1))
  if(sum(clone_cnv!=1) < cnv_thre){
    next
  }
  copy_ratio_pure <- rep(1, length(clone_cnv))
  f <- clone_cnv > 1
  copy_ratio_pure[f] <- copy_ratio_pure_gain[f]
  f <- clone_cnv < 1
  copy_ratio_pure[f] <- copy_ratio_pure_loss[f]
  
  clone_cell <- cell_barcodes[clone == k]
  tumor_content_estimate <- rep(0, length(clone_cell))
  for(i in 1:length(clone_cell)){
    f <- clone_cnv != 1
    df_ <- data.frame(x = copy_ratio_pure[f] - 1, 
                      y = copy_ratio_seg[clone_cell[i],f] - 1)
    re <- lm(data = df_, formula = y ~ x + 0)
    df[clone_cell[i], "p_estimate"] <- re$coefficients
    df[clone_cell[i], "p_estimate_confidence"] <- summary(re)$adj.r.squared
  }
}

## visualize
obj <- AddMetaData(obj, df$p_estimate, col.name="spacnv")
SpatialFeaturePlot(obj, features = "spacnv")