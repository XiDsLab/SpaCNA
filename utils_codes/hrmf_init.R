## Chromosome information of simulated data cnv=3
cnv_setting_cnv3 <- data.frame(chr=paste0("chr", c(1:22)), cnv=1)
# cnv_setting_cnv3$cnv[c(3,4,22)] <- 0.5 #3
# cnv_setting_cnv3$cnv[c(5,6,7,20,21)] <- 0.5 #5
# cnv_setting_cnv3$cnv[c(9,10,11,15,16)] <- 1.5 # 5
# cnv_setting_cnv3$cnv[c(12,13,14)] <- 1.5 #3

cnv_setting_cnv3$cnv[c(1,3,5)] <- 0.5
cnv_setting_cnv3$cnv[c(7,9,11,13)] <- 1.5

## simulate the expression matrix
sim_mat_new<- function(expr,gene_means,num_cells, common_dispersion,dropout_prob='au') {
  ##if dropout_prob='au' get the parameters of 0 else use the given number
  ngenes <- length(gene_means)
  sim_matrix <- matrix(rep(0,ngenes*num_cells), nrow=ngenes)
  
  ## Simulation performed on each gene
  if (dropout_prob=="au"){
    m <- rowMeans(expr)
    p0 <- apply(expr, 1, function(x) { sum(x==0) })/ncol(expr)
    para_0 <- smooth.spline(log(m[m>0]), p0[m>0])
    
    sim_expr_vals <- function(gene_idx) {
      gene_mean <- gene_means[gene_idx]
      val <- 0
      if (gene_mean > 0) {
        val <- rnbinom(n = 1,
                       mu = gene_mean,
                       size = 1 / common_dispersion)
        if (val > 0) {
          dropout_prob <- predict(para_0, log(val))$y[1]
          if (runif(1) <= dropout_prob) {
            val <- 0
          }
        }
      }
      return(val)
    }
  } else{
    sim_expr_vals <- function(gene_idx) {
      gene_mean <- gene_means[gene_idx]
      val <- 0
      if (gene_mean > 0) {
        val <- rnbinom(n = 1,
                       mu = gene_mean,
                       size = 1 / common_dispersion)
        if (val > 0) {
          if (runif(1) <= dropout_prob) {
            val <- 0
          }
        }
      }
      return(val)
    }
    
  }
  
  ## get the simulated matrix
  for (i in 1:num_cells) {
    sim_vals <- sapply(1:ngenes, FUN=sim_expr_vals)
    sim_matrix[,i] <- sim_vals
  }
  
  return(sim_matrix)
}

get_sim_expr_new <- function(expr, 
                             normal_idx, 
                             gene_pos, 
                             gene_normalize=T,
                             cnv_num=3, 
                             num_cells=100,
                             tumor_perc=1,
                             common_dispersion=0.1,
                             dropout_prob='au'
                         ) {
  if(cnv_num==3){
    cnv_setting=cnv_setting_cnv3
  }else {
    cnv_setting=cnv_setting_cnv5
  }
  ####intersect
  intersection<-intersect(rownames(expr),gene_pos$symbol)
  gene_pos<-gene_pos[gene_pos$symbol %in% intersection,]
  expr<-expr[gene_pos$symbol,]
  ## normalize 
  if(gene_normalize){
    cs <- colSums(expr)
    expr <- sweep(expr, STATS=cs, MARGIN=2, FUN="/")*median(cs)
  }

  expr_normal <- expr[,normal_idx]
  gene_means <- rowMeans(expr_normal)
  gene_means[gene_means==0] <- 1e-3 
  names(gene_means) <- gene_pos$symbol       
  
  ## simulate normal cells:
  print("Simulating normal cells...")
  sim_normal_matrix <- sim_mat_new(expr,
                               gene_means,
                               num_cells,
                               common_dispersion=common_dispersion,
                               dropout_prob=dropout_prob)
  #print(nrow(sim_normal_matrix))
  #print(length(gene_pos$symbol))
  rownames(sim_normal_matrix) <- gene_pos$symbol
  colnames(sim_normal_matrix) <- paste0("normal", seq_len(num_cells))
  ## simulate tumor cells
  print("Simulating tumor cells...")
  tumor_gene_means <- gene_means
  #print(nrow(cnv_setting))
  for (i in 1:nrow(cnv_setting)) {
    chr <- cnv_setting$chr[i]
    cnv <- cnv_setting$cnv[i]
    if (cnv != 1) {
      gene_idx <- which(gene_pos$chr==chr)
      tumor_gene_means[gene_idx] <-  tumor_gene_means[gene_idx]*cnv
    }
  }
  sim_tumor_matrix <- sim_mat_new(expr, 
                        tumor_perc*tumor_gene_means,
                        num_cells,
                        common_dispersion=1/tumor_perc*common_dispersion,
                        dropout_prob=dropout_prob)

  
  if (tumor_perc!=1){
    #print(1/(1-tumor_perc)*common_dispersion)
    sim_tumor_immu_matrix <- sim_mat_new(expr, 
                                (1-tumor_perc)*gene_means,
                                num_cells,
                                common_dispersion=1/(1-tumor_perc)*common_dispersion,
                                dropout_prob=dropout_prob)
    sim_tumor_matrix<-sim_tumor_matrix+sim_tumor_immu_matrix

  }
  rownames(sim_tumor_matrix) <- gene_pos$symbol
  colnames(sim_tumor_matrix) <- paste0("tumor", seq_len(num_cells))  
  sim_expr <- cbind(sim_normal_matrix, sim_tumor_matrix) 
  cell_type <- c(rep("normal", num_cells), rep("tumor", num_cells))
  
  sim_info<- list(sim_expr=sim_expr,
                  cell_type=cell_type,
                  cnv_setting=cnv_setting,
                  tumor_perc=tumor_perc)
  return(sim_info)
}

estimate_hmrf_parameter <- function(count_RNA, 
                                    normal_cells,
                                    gene_pos,
                                    bin,
                                    gene_normalize=F,
                                    cnv_num=3,
                                    num_cells=200,
                                    tumor_perc=1,
                                    common_dispersion=0.1,
                                    dropout_prob="au",
                                    dlm_dV=0.1, 
                                    dlm_dW=0.001,
                                    alpha=c(0.5, 0.25, 0.25),
                                    gene_thre=2){
  to_use_info <- get_sim_expr_new(count_RNA[,normal_cells], 
                                  normal_idx=normal_cells, 
                                  gene_pos=gene_pos,
                                  gene_normalize=gene_normalize,
                                  cnv_num=cnv_num,
                                  num_cells=num_cells,
                                  tumor_perc=tumor_perc,
                                  common_dispersion=common_dispersion,
                                  dropout_prob=dropout_prob)
  
  count_sim_RNA <- edgeR::cpm(to_use_info$sim_expr)
  count_sim_smooth <- chr_smooth(count_sim_RNA, log=T, 
                                 dlm_dV=dlm_dV, dlm_dW=dlm_dW)
  count_sim_smooth <- edgeR::cpm(count_sim_smooth)
  neigh_expr <- get_expr_neigh_list(count_sim_smooth, pca=20, k=6,
                                    thre=0.5)
  count_sim_smooth <- spatial_smooth(count_sim_smooth,
                                     neigh_expr,
                                     neigh_expr,
                                     alpha=alpha)
  count_sim_norm <- get_copy_ratio(count_sim_smooth, 
                                   to_use_info$cell_type=="normal",
                                   thre="none")
  
  count_sim_norm <- generate_bin_count(count_sim_norm, 
                                       gene_pos, 
                                       bin, 
                                       func=mean)
  gene_num <- rep(0, nrow(bin))
  for(i in 1:nrow(bin)){
    f <- gene_pos$chr==bin$chr[i] & gene_pos$start>=bin$start[i] & gene_pos$end<=bin$end[i]
    gene_num[i] <- sum(f)
  }
  count_sim_norm_2 <- count_sim_norm[gene_num>gene_thre, ]
  state_ratio <- c(0.5,1,1.5)
  mus <- rep(0, 3)
  sigmas <- rep(0, 3)
  for(i in 1:3){
    cnv_chr <- to_use_info$cnv_setting$chr[to_use_info$cnv_setting$cnv==state_ratio[i]]
    f1 <- str2bin(rownames(count_sim_norm_2))$chr %in% cnv_chr
    mus[i] <- mean(count_sim_norm_2[f1, to_use_info$cell_type=="tumor"])
    sigmas[i] <- sd(count_sim_norm_2[f1, to_use_info$cell_type=="tumor"])
  }
  return(list(mus, sigmas))
}







