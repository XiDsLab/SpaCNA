generate_bin_count <- function(count_RNA, gene_pos, bin, func=mean){
  count_RNA <- count_RNA[gene_pos$symbol,]
  idx <- duplicated(rownames(count_RNA))
  count_RNA <- count_RNA[!idx,]
  gene_pos <- gene_pos[!idx,]
  
  count_bin <- matrix(nrow=nrow(bin), ncol=ncol(count_RNA))
  for(i in 1:nrow(bin)){
    # print(i)
    f <- gene_pos$chr==bin$chr[i] & gene_pos$start>=bin$start[i] & gene_pos$end<=bin$end[i]
    f <- gene_pos$symbol[f]
    if(length(f)==0){
      count_bin[i,] <- 0
    }else if(length(f)==1){
      count_bin[i,] <- count_RNA[f,]
    }else{
      count_bin[i,] <- apply(count_RNA[f,], 2, func)
    }
  }
  colnames(count_bin) <- colnames(count_RNA)
  rownames(count_bin) <- bin2str(bin)
  return(count_bin)
  # saveRDS(count_bin, paste0(sample_dir,"count_bin.rds"))
}

preprocess <- function(count, bin, normal_cells="none", 
                       thre_bin=0.8, thre_cell=1, 
                       thre_map=500000){
  # count: bin*cell
  count_filter <- count
  temp <- colSums(count_filter<=0)
  f1 <- temp<thre_cell*nrow(count_filter)
  print(paste("number of cell before filter:", length(f1)))
  print(paste("number of cell after filter:", sum(f1)))
  count_filter <- count_filter[,f1]
  
  if(identical(normal_cells, "none")){
    temp <- rowSums(count_filter<=0)
    f2 <- temp<thre_bin*ncol(count_filter) & bin$map>=thre_map
  }else{
    temp <- rowSums(count_filter[,normal_cells]<=0)
    f2 <- temp<thre_bin*ncol(count_filter[,normal_cells]) & bin$map>=thre_map
  }

  print(paste("number of bin before filter:", length(f2)))
  print(paste("number of bin after filter:", sum(f2)))
  count_filter <- count_filter[f2,]
  
  q <- quantile(count_filter, probs=c(1:100)*0.01)
  thre <- q[99]
  count_filter[count_filter>thre] <- thre
  return(count_filter)
}


get_copy_ratio <- function(count, normal_cell, thre="none"){
  baseline <- apply(count[,normal_cell], 1, median)
  norm_count <- count/baseline
  if(thre == "none"){
    q <- quantile(norm_count, probs=c(1:100)*0.01)
    print(paste0("99th percentile: ",q[99]))
    thre <- q[99]
  }
  norm_count[norm_count>thre] <- thre
  return(norm_count)
}


get_expr_neigh_list<-function(expr,pca=20,k=10,thre=0.4){
  
  neigh_list<-list()
  pca_result <- prcomp(t(expr), scale = TRUE,center = TRUE,rank = pca)
  feature_pca <- pca_result$x
  feature_cor<- cor(t(feature_pca),t(feature_pca))
  feature_neigh<-feature_cor>=thre
  
  knn<- FindNeighbors(as.matrix(feature_pca), 
                      k.param=k, 
                      return.neighbor=TRUE)
  for(i in 1:nrow(feature_pca)){
    idx1 <- knn@nn.idx[i,-1] 
    idx2<- feature_neigh[i,idx1] 
    neigh_final<-idx1[idx2]
    neigh_list[[i]]<- neigh_final
  }
  return(neigh_list)
  
}


chr_smooth <- function(count, log=T, dlm_dV=0.1, dlm_dW=0.02){
  q <- quantile(count, probs=c(1:100)*0.01)
  thre <- q[99]
  count[count>thre] <- thre
  
  count <- edgeR::cpm(count)
  
  if(log){
    count <- log(count+1)
  }else{
    count <- count
  }
  dlm_sm <- function(y){
    model <- dlm::dlmModPoly(order=1, dV=dlm_dV, dW=dlm_dW)
    x <- dlm::dlmSmooth(y, model)$s
    x <- x[2:length(x)]
    return(x)
  }
  temp <- apply(count, 2, dlm_sm)
  if(log){
    count <- exp(temp)
  }else{
    count <- temp
  }
  return(count)
}


spatial_smooth <- function(count, neigh_expr, neigh_cell,
                           alpha=c(0.5, 0.25, 0.25)){
  count_smooth_spa<-count
  for(i in 1:ncol(count)){
  idx_expr<-neigh_expr[[i]]
  if(length(idx_expr)>=2){
    expr_nei<-rowMeans(count[,idx_expr])
  }else if(length(idx_expr)==1){
    expr_nei<-count[,idx_expr]
  }else{
    expr_nei<-count[,i]
  }
  
  idx_cell<-neigh_cell[[i]]
  if(length(idx_cell)>=2){
    cell_nei<-rowMeans(count[,idx_cell])
  }else if(length(idx_cell)==1){
    cell_nei<-count[,idx_cell]
  }else{
    cell_nei<-count[,i]
  }
  count_smooth_spa[,i]<-alpha[1]*count[,i]+ 
                        alpha[2]* expr_nei+ 
                        alpha[3]*cell_nei
  }
  return(count_smooth_spa)
}


