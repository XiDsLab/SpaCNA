library(Seurat)
library(biomaRt)
library(ComplexHeatmap)
library(parallelDist)
library(irlba)
library(edgeR)
library(ggplot2)
library(rootSolve)
library(patchwork)
library(glmnet)
#library(reticulate)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

merge_bkp <- function(bkp1,bkp2,min_interval=3){
  # bkp2 all reserved
  for(i in bkp1){
    if(min(abs(bkp2-i))>min_interval){
      bkp2[length(bkp2)+1] <- i
    }
  }
  return(bkp2)
}

str2bin <- function(bin_name){
  temp <- strsplit(bin_name, split="_")
  temp2 <- unlist(temp)
  dim(temp2) <- c(3,length(temp2)/3)
  bins <- data.frame(chr=temp2[1,],start=as.numeric(temp2[2,]),
                     end=as.numeric(temp2[3,]))
  bins$n <- c(1:nrow(bins))
  return(bins)
}

bin2str <- function(bin){
  return(paste0(bin$chr,"_",bin$start,"_",bin$end))
}

bin2chrbkp <- function(bins){
  chr_bkp <- c(0)
  temp <- unique(bins$chr)
  for(i in 1:length(temp)){
    chr_bkp[i+1] <- max(bins$n[bins$chr==temp[i]])
  }
  return(chr_bkp)
}

del_chr <- function(count, chr_list=c("chrX","chrY")){
  bin <- str2bin(colnames(count))
  return(count[,!bin$chr %in% chr_list])
}


get_gene_pos <- function(gene_list){
  ensembl <- useMart("ensembl",host="https://www.ensembl.org")
  #ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                   #path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  
  ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
  gene_pos <- getBM(attributes=c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
                    filters=c('hgnc_symbol'),
                    values=list(gene_list),
                    mart=ensembl)
  gene_pos <- gene_pos[gene_pos$chromosome_name %in% c(1:22,"X","Y"),]
  colnames(gene_pos) <- c("symbol","chr","start","end")
  gene_pos$chr <- paste0("chr",gene_pos$chr)
  return(gene_pos)
}

get_spot_location <- function(seurat_obj, type="10X"){
  if(type=="10X"){
    spot_location <- seurat_obj@images[["slice1"]]@coordinates
    spot_location <- spot_location[,c("imagerow","imagecol")]
    names(spot_location) <- c("x","y")
    return(spot_location)
  }else{
    spot_location <- seurat_obj@images[["image"]]@coordinates
    spot_location <- spot_location[,c("x","y")]
    return(spot_location)
  }
}

view_distribution <- function(r, binwidth=0.020){
  df <- data.frame(r=r)
  p <- ggplot(df, aes(x=r))+geom_histogram(binwidth=binwidth)
  print(p)
}

plot_mix_comps <- function(x, mu, sigma, lambda){
  return(lambda*dnorm(x, mu, sigma))
}

plot_aggre_comps <- function(x, mu, sigma, lambda){
  y <- 0
  for(i in 1:length(mu)){
    y <- y+lambda[i]*dnorm(x, mu[i], sigma[i])
  }
  return(y)
}

plot_GMM <- function(r, binwidth=0.02, mu, sigma, lambda, aggre=T){
  df <- data.frame(r=r)
  p <- ggplot(df)+geom_histogram(aes(x=r, y=..density..), binwidth=0.01)
  if(aggre){
    p <- p+stat_function(geom="line", fun=plot_aggre_comps,
                         args=list(mu, sigma, lambda))
  }
  if(!aggre){
    for(i in 1:length(mu)){
      p <- p+stat_function(geom="line", fun=plot_mix_comps,
                           args=list(mu[i], sigma[i], lambda[i]))
    }
  }
  print(p)
}

denoise<-function(mat,thre=0.95){
  mat2<-mat
  n_bin<-ncol(mat)
  for (i in c(1:nrow(mat))){
    temp<-mat[i,]
    most<-Mode(temp)
    if(sum(temp==most)/n_bin >=thre){
      mat2[i,]<-most
    }
  }
  return(mat2)
}


generate_bin_cnv<- function(gene_pos, bin){
  bin_new<-bin
  for(i in 1:nrow(bin)){
    #print(i)
    f <- gene_pos$chr==bin$chr[i] & gene_pos$start>=bin$start[i] & gene_pos$end<=bin$end[i]
    f <- gene_pos$cnv[f]
    bin_new$cnv[i]<-Mode(f)
  }
  return(na.omit(bin_new))
}


get_neigh_list<-function(feature_pca,spot_location,k=7,thre=0.2,dist_thre=500){
  
  
  feature_pca<-feature_pca[rownames(spot_location),]
  neigh_list<-list()
  #adj <- matrix(0, nrow=nrow(spot_location), ncol=nrow(spot_location))
  #feature_pca<-feature_pca[rownames(spot_location),]
  feature_cor<- cor(t(feature_pca),t(feature_pca))
  feature_neigh<-feature_cor>=thre
  
  knn<- FindNeighbors(as.matrix(spot_location), 
                      k.param=k, 
                      annoy.metric="manhattan", 
                      return.neighbor=TRUE)
  for(i in 1:nrow(spot_location)){
    idx1 <- knn@nn.idx[i,-1] #######spatial neighbour
    dist1<-knn@nn.dist[i,-1]
    idx1 <- idx1[dist1<=dist_thre] ####remove spot farther
    idx2<- feature_neigh[i,idx1] ####feature neighbour
    neigh_final<-idx1[idx2]
    neigh_list[[i]]<- neigh_final
  }
  names(neigh_list)<-rownames(spot_location)
  return(neigh_list)
  
}

get_neigh_list_simple <- function(spot_location, k = 6, dist_thre = 500){
  neigh_list<-list()
  knn<- FindNeighbors(as.matrix(spot_location), 
                      k.param=k, 
                      annoy.metric="manhattan", 
                      return.neighbor=TRUE)
  for(i in 1:nrow(spot_location)){
    idx1 <- knn@nn.idx[i,-1] 
    dist1 <- knn@nn.dist[i,-1]
    idx1 <- idx1[dist1<=dist_thre] 
    neigh_list[[i]]<- idx1
  }
  names(neigh_list)<-rownames(spot_location)
  return(neigh_list)
}

calculate_moran_I <- function(y, W) {
  y_bar <- mean(y)
  
  N <- length(y)
  W_sum <- sum(W)

  numerator <- sum(W * outer(y - y_bar, y - y_bar))
  denominator <- sum((y - y_bar)^2)
  I <- (N / W_sum) * (numerator / denominator)
  return(I)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

get_final_neigh<-function(image_dir,spot_location,k=7,thre=0.2,dist_thre='unknown'){
  image_feature<- read.table(paste0(image_dir,'/resnet50/features_pca.txt'))
  image_name<- read.table(paste0(image_dir,'/resnet50/feature_name.txt'))
  rownames(image_feature)<-image_name$V1
  if(is.numeric(dist_thre)){
    dist_thre<-dist_thre
  }else{
    knn2 <- FindNeighbors(as.matrix(spot_location), 
                          k.param=7, 
                          annoy.metric="manhattan", 
                          return.neighbor=TRUE)
    ###check the threshold
    dist_thre<-Mode(knn2@nn.dist[c(1:nrow(spot_location)),6])
  }
  neigh_list<- get_neigh_list(image_feature,spot_location,k,thre,dist_thre)
  return(neigh_list)
  
}
save_image_info<-function(obj,image_dir,type="10X"){
  if(type=="10X"){
    scalefactors<-obj@images$slice1@scale.factors$hires
    coordinate<-obj@images$slice1@coordinates[4:5]*scalefactors
    spot_name<-rownames(coordinate)
  }else{
    scalefactors<-obj@images$image@scale.factors$hires
    coordinate<-obj@images$image@coordinates[4:5]*scalefactors
    spot_name<-rownames(coordinate)
  }
  write.table (coordinate, file =paste0(image_dir, "exp_location.txt"), sep =" ", row.names =FALSE, col.names =FALSE, quote =TRUE)
  write.table (spot_name, file =paste0(image_dir, "spot.txt"), sep =" ", row.names =FALSE, col.names =FALSE, quote =TRUE)
  knn1<- FindNeighbors(as.matrix(coordinate), 
                       k.param=2, 
                       annoy.metric="manhattan", 
                       return.neighbor=TRUE)
  r<-1/2*Mode(knn1@nn.dist[c(1:nrow(coordinate)),-1])
  
  return(list(coordinate=coordinate,spot_name=spot_name,r=r))
}

get_image_feature<-function(py_env,py_dir,image_dir,r=15){
  use_condaenv(py_env)
  source_python(py_dir)
  get_pca_feature(image_dir,r=r)
  image_feature<- read.table(paste0(image_dir,'/resnet50/features_pca.txt'))
  image_name<- read.table(paste0(image_dir,'/resnet50/feature_name.txt'))
  rownames(image_feature)<-image_name$V1
  return(image_feature)
  
}

