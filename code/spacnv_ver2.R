source("./code_ver2/utils.R")
source("./code_ver2/downstream.R")
source("./code_ver2/normalize.R")
source("./code_ver2/hrmf_em.R")
source("./code_ver2/hrmf_init.R")

sample_dir <- "./data/CRC_2/"
plot_dir <- "./CRC_2_analysis/"
obj <- readRDS("./data/CRC_2/seurat_object.rds")

## inputs
count_RNA <- as.matrix(obj@assays$Spatial@counts)
count_RNA <- count_RNA[rowMeans(count_RNA)>0.05, ]
print(dim(count_RNA))
bin <- readRDS("./data/bin_K=1e+06_hg38.RDS")
cell_barcodes <- colnames(count_RNA)
seurat_cluster <- as.numeric(obj$seurat_clusters)-1
cluster <- as.character(seurat_cluster)
spot_location <- get_spot_location(obj)
normal_cells <- cell_barcodes[! cluster %in% c(7, 12)]

# gene position should be ordered 
gene_pos <- readRDS(paste0(sample_dir, "gene_pos.rds"))
gene_pos <- gene_pos[gene_pos$symbol %in% rownames(count_RNA), ]
gene_pos$n <- 1:nrow(gene_pos)

# count gene number in bins
gene_num <- rep(0, nrow(bin))
for(i in 1:nrow(bin)){
  f <- gene_pos$chr==bin$chr[i] & gene_pos$start>=bin$start[i] & gene_pos$end<=bin$end[i]
  gene_num[i] <- sum(f)
}

## gene smooth
count_RNA <- count_RNA[gene_pos$symbol,]

count_RNA_smooth <- chr_smooth(count_RNA, log=T, dlm_dV=0.1, dlm_dW=0.001)
count_RNA_smooth <- edgeR::cpm(count_RNA_smooth)
neigh_cell <- get_neigh_list_simple(spot_location)
count_RNA_spa <- spatial_smooth(count_RNA_smooth,
                                neigh_cell,
                                neigh_cell,
                                alpha=c(0.5, 0.25, 0.25))
count_RNA_norm <- get_copy_ratio(count_RNA_spa, normal_cells)

# f <- sample(1:length(cluster), 1000)
# plot_heatmap(t(count_RNA_norm)[f,], 
#              bin=gene_pos,
#              cell_cluster=cluster[f],
#              cell_annotation=cluster[f],
#              output_dir=plot_dir,
#              output_name="copy_ratio_gene.png")
# obj <- AddMetaData(obj, colMeans(abs(count_RNA_norm-1)), col.name="cnv_score")
# SpatialFeaturePlot(obj, features="cnv_score")

count_norm <- generate_bin_count(count_RNA_norm, gene_pos, bin, func=mean)
# filter out bins with few genes
f <- gene_num >= 2
count_norm <- count_norm[f,]
saveRDS(count_norm, paste0(sample_dir, "count_norm.rds"))

# plot
f <- sample(1:length(cluster), 1000)
plot_heatmap(t(count_norm)[f,], 
             cell_cluster=cluster[f],
             cell_annotation=cluster[f],
             output_dir=plot_dir,
             output_name="copy_ratio_bin.png")
# obj <- AddMetaData(obj, colMeans(abs(count_norm-1)), col.name="cnv_score")
# SpatialFeaturePlot(obj, features = "cnv_score")

## hmrf para init
parameter_init <- estimate_hmrf_parameter(count_RNA=count_RNA, 
                                          normal_cells=normal_cells,
                                          gene_pos=gene_pos,
                                          bin=bin,
                                          gene_normalize=F,
                                          cnv_num=3,
                                          num_cells=200,
                                          tumor_perc=1,
                                          common_dispersion=0.1,
                                          dropout_prob="au",
                                          dlm_dV=0.1, 
                                          dlm_dW=0.001,
                                          alpha=c(0.5, 0.25, 0.25),
                                          gene_thre=2)
saveRDS(parameter_init, paste0(sample_dir, "hrmf_para.rds"))

## cn state init
count_norm <- readRDS(paste0(sample_dir, "count_norm.rds"))
parameter_init <- readRDS(paste0(sample_dir, "hrmf_para.rds"))
mus <- parameter_init[[1]]
sigmas <- parameter_init[[2]]
mus_df <- data.frame(row.names=cell_barcodes,
                     p=rep(1, length(cell_barcodes)),
                     mu1=mus[1], mu2=mus[2], mu3=mus[3])
sigs_df <- data.frame(row.names=cell_barcodes,
                      p=rep(1, length(cell_barcodes)),
                      mu1=sigmas[1], mu2=sigmas[2], mu3=sigmas[3])

cnv_state <- get_cnv_init_perce(count_norm, mus_df, sigs_df)
state_ratio <- c(0.5, 1, 1.5)
cnv_state_init <- state2ratio(cnv_state, state_ratio)
f <- sample(1:length(cluster), 1000)
plot_heatmap(t(cnv_state_init)[f,],
             cell_annotation=cluster[f],
             output_dir=plot_dir,
             output_name="cns_init_ver2.png")

## segment
baseline <- rowMeans(count_norm[,normal_cells])
bk_bic <- calculate_CNV(t(count_norm), 
                        baseline, 
                        lambda=0.0075,
                        sample_dir=sample_dir)
saveRDS(bk_bic, paste0(sample_dir, "bk_bic.rds"))

## hmrf
cns <- get_state_icm_seg_perc(count_norm, 
                              bk_bic, 
                              neigh_list=neigh_cell, 
                              mus_df=mus_df, 
                              sigmas_df=sigs_df, 
                              max_iter=5,
                              beta_fixed=1)
cns_ <- cns
cns <- state2ratio(cns_, state_ratio)
saveRDS(cns, paste0(sample_dir, "cns.rds"))

cns <- denoise(cns, thre=0.9)
f <- sample(1:length(cluster), 1000)
plot_heatmap(t(cns)[f,],
             cell_annotation=cluster[f],
             output_dir=plot_dir,
             output_name="cns_ver2.png")
obj <- AddMetaData(obj, colSums(cns!=1), col.name="cnv_score")
SpatialFeaturePlot(obj, features = "cnv_score")

#### tumor content ####
## cell2location 
deconv_re <- readRDS(paste0(sample_dir, "cell2location_re.rds"))
obj <- AddMetaData(obj, deconv_re, col.name="cell2location")
SpatialFeaturePlot(obj, features = "cell2location")

## spacnv
cns <- readRDS(paste0(sample_dir, "cns.rds"))
count_norm <- readRDS(paste0(sample_dir, "count_norm.rds"))
bk_bic <- readRDS(paste0(sample_dir, "bk_bic.rds"))
cell_barcodes <- colnames(count_norm)
K <- 7
cnv_thre <- 10

temp <- count_norm
d <- parallelDist::parallelDist(t(temp))
hc_re <- hclust(d, method="ward.D2")
clone <- cutree(hc_re, k=K)
sum(clone == 5)

obj <- AddMetaData(obj, clone, col.name="clone")
SpatialDimPlot(obj, group="clone")

plot_heatmap(t(cns),
             cell_cluster=clone,
             cell_annotation=clone,
             output_dir=plot_dir,
             output_name="cns_clone_annot_ver2.png")

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

tumor_content <- data.frame(row.names = colnames(cns),
                            p <- rep(0, ncol(cns)))
df <- tumor_content
df$p_estimate <- 0
df$p_estimate_confidence <- 0

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
  }else if(sum(x == 1.5)/length(x) > thre){
    return(1.5)
  }else{
    return(1)
  }
}

for(k in c(1:K)){
  # print(sum(clone==k))
  clone_cnv <- apply(cns_seg[clone==k,], 2, get_cnv_clone)
  print(sum(clone_cnv!=1))
  if(sum(clone_cnv!=1) < cnv_thre){
    next
  }
  # copy_ratio_pure <- rep(1, length(clone_cnv))
  # for(i in 1:length(clone_cnv)){
  #   if(clone_cnv[i] == 0.5){
  #     copy_ratio_pure[i] <- mean(sort(copy_ratio_seg[clone==k, i])[1:5])
  #   }
  #   if(clone_cnv[i] == 1.5){
  #     copy_ratio_pure[i] <- mean(sort(copy_ratio_seg[clone==k, i], 
  #                                     decreasing = T)[1:5])
  #   }
  # }
  copy_ratio_pure <- rep(1, length(clone_cnv))
  f <- clone_cnv > 1
  copy_ratio_pure[f] <- copy_ratio_pure_gain[f]
  f <- clone_cnv < 1
  copy_ratio_pure[f] <- copy_ratio_pure_loss[f]
  
  clone_cell <- cell_barcodes[clone == k]
  tumor_content_estimate <- rep(0, length(clone_cell))
  for(i in 1:length(clone_cell)){
    # f <- cns_seg[clone_cell[i],] != 1 & clone_cnv != 1
    f <- clone_cnv != 1
    df_ <- data.frame(x = copy_ratio_pure[f] - 1, 
                      y = copy_ratio_seg[clone_cell[i],f] - 1)
    # plot(df$y, df$x)
    re <- lm(data = df_, formula = y ~ x + 0)
    df[clone_cell[i], "p_estimate"] <- re$coefficients
    df[clone_cell[i], "p_estimate_confidence"] <- summary(re)$adj.r.squared
  }
}

ggplot(df, aes(x=p, y=p_estimate)) + geom_point()
obj <- AddMetaData(obj, df$p_estimate, col.name="spacnv")
SpatialFeaturePlot(obj, features = "spacnv")

## spacet
spacet_re <- readRDS(paste0(sample_dir, "spacet_re.rds"))
obj <- AddMetaData(obj, as.vector(spacet_re[1,]), col.name="spacet")
SpatialFeaturePlot(obj, features="spacet")

## infercnv
infercnv_re <- readRDS(paste0(sample_dir, "infercnv_re"))
infercnv_re <- infercnv_re@expr.data
obj <- AddMetaData(obj, colSums(infercnv_re!=2), col.name="infercnv")
SpatialFeaturePlot(obj, features="infercnv")

## plot
png(paste0(plot_dir, "tumor_content.png"), width=2000, height=2000)
SpatialFeaturePlot(obj, features = c("spacnv", "infercnv", 
                                     "cell2location", "spacet"),
                   ncol=2) & 
  theme(legend.title=element_text(size=50))
dev.off()

saveRDS(obj, paste0(sample_dir, "seurat_object.rds"))

##
df <- data.frame(x=obj$spacnv, y=obj$cell2location)
png(paste0(plot_dir, "spacnv_vs_cell2location.png"), width=500, height=500)
ggplot(df, aes(x=x, y=y)) + geom_point() + theme_bw()
dev.off()

df <- data.frame(x=obj$spacnv, y=obj$spacet)
png(paste0(plot_dir, "spacnv_vs_spacet.png"), width=500, height=500)
ggplot(df, aes(x=x, y=y)) + geom_point() + theme_bw()
dev.off()