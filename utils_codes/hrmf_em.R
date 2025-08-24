emission_probs <- function(x, mus, sigmas){
  p <- c(1:length(mus))
  for(i in 1:length(p)){
    p[i] <- dnorm(x,mus[i],sigmas[i])
  }
  return(p)
}


cnv_state_init <- function(data, mus, sigmas){
  emission_probs <- sapply(data, emission_probs, 
                           mus=mus, sigmas=sigmas, simplify=TRUE)
  R <- apply(emission_probs, 2, which.max)
  return(R)
}


emission_probs_matrix <- function(norm_data,mus_df,sigmas_df){
  # norm_data: bin*cell
  n_state <- ncol(mus_df) - 1
  emission_prob_mat <- array(NA, dim = c(nrow(norm_data), ncol(norm_data), 
                                         n_state))
  for(j in 1:ncol(norm_data)){
    for(k in 1:n_state){
      emission_prob_mat[, j, k] <- dnorm(norm_data[,j], mus_df[j,k + 1], sigmas_df[j,k + 1])
    }
  }
  return(emission_prob_mat)
}

get_cnv_init_perce<-function(norm_data,mus_df,sigmas_df){
  ## Align cell names
  mus_df<-mus_df[colnames(norm_data),]
  sigmas_df<-sigmas_df[colnames(norm_data),]
  emission_prob_mat <- emission_probs_matrix(norm_data, mus_df, sigmas_df)
  cnv_state_init <- apply(emission_prob_mat, c(1,2), which.max)
  rownames(cnv_state_init) <- rownames(norm_data)
  colnames(cnv_state_init) <- colnames(norm_data)
  return(cnv_state_init)
}


get_beta <- function(cnv_state, neigh_list, beta_range=c(0.1,3)){
  state <- unique(cnv_state)
  state_num <- length(state)
  if(state_num==1){
    return(-1)
  }
  # print(state_num)
  cnv_neigh <- matrix(1, length(cnv_state), state_num+1)
  ## column j: the number of neighbors labeled jth state
  ## column state_num+1: the number of neighbors labeled cnv_i
  for (i in (1:length(cnv_state))){
    neigh_ <- neigh_list[[i]]
    j <- 1
    for (k in state){
      cnv_neigh[i,j] <- sum(cnv_state[neigh_]==k)
      j <- j+1
    }
    cnv_neigh[i,state_num+1] <- sum(cnv_state[neigh_]==cnv_state[i])
  }
  ## define the beta-function
  fun <- function(beta){
    if(state_num>1){
      tmp1 <- rowSums(cnv_neigh[,1:state_num]*exp(beta*cnv_neigh[,1:state_num]))
      tmp2 <- rowSums(exp(beta*cnv_neigh[,1:state_num]))
    }else{
      tmp1 <- cnv_neigh[,1:state_num]*exp(beta*cnv_neigh[,1:state_num])
      tmp2 <- exp(beta*cnv_neigh[,1:state_num])
    }
    return(sum(cnv_neigh[,state_num+1])-sum(tmp1/tmp2))
  }
  ## calculate the zero point
  beta <- uniroot.all(Vectorize(fun), beta_range)
  if(length(beta)==0){
    return(-1)
  }else if(length(beta)>1){
    return(beta[1])
  }else{
    return(beta)
  }
}


get_beta_seg <- function(cnv_state, neigh_list, beta_range=c(0.1,3)){
  state <- unique(c(cnv_state))
  state_num <- length(state)
  if(state_num==1){
    return(-1)
  }
  # print(state_num)
  cnv_neigh <- array(1, c(nrow(cnv_state),ncol(cnv_state) ,state_num+1))
  ## column j: the number of neighbors labeled jth state
  ## column state_num+1: the number of neighbors labeled cnv_i
  for (i in (1:nrow(cnv_state))){
    for (j in (1:ncol(cnv_state))){
      neigh_ <- neigh_list[[j]]
      d<-1
      for (k in state){
        #print(neigh_)
        cnv_neigh[i,j,d] <- sum(cnv_state[i,neigh_]==k)
        d <- d+1
      }
      cnv_neigh[i,j,state_num+1] <- sum(cnv_state[i,neigh_]==cnv_state[i,j])
    }
  }
  
  ## define the beta-function
  fun <- function(beta){
    if(state_num>1){
      tmp1 <- apply(cnv_neigh[,,1:state_num]*exp(beta*cnv_neigh[,,1:state_num]),c(1:2),sum)
      tmp2 <- apply(exp(beta*cnv_neigh[,,1:state_num]),c(1:2),sum)
    }else{
      tmp1 <- cnv_neigh[,,1:state_num]*exp(beta*cnv_neigh[,,1:state_num])
      tmp2 <- exp(beta*cnv_neigh[,,1:state_num])
    }
    return(sum(sum(cnv_neigh[,,state_num+1]))-sum(sum(tmp1/tmp2)))
  }
  ## calculate the zero point
  beta <- uniroot.all(Vectorize(fun), beta_range)
  if(length(beta)==0){
    return(-1)
  }else if(length(beta)>1){
    return(beta[1])
  }else{
    return(beta)
  }
}

get_state_icm_seg_perc<-function(expr, breaks_union,neigh_list, mus_df, sigmas_df, beta_fixed='au',beta_default=1, max_iter=5){
  t1 <- Sys.time()
  n_spot <- ncol(expr)
  n_state <- ncol(mus_df)-1
  
  breaks_union<-sort(breaks_union)
  breaks_union[length(breaks_union)]<-breaks_union[length(breaks_union)]+1
  n_break<-length(breaks_union)
  emission_prob_mat <- emission_probs_matrix(expr, mus_df, sigmas_df)
  cnv_state_init <- apply(emission_prob_mat, c(1,2), which.max)
  if (any(is.na(cnv_state_init))) {
    print("cnv_state_init contains NA values!")
}

  mus<-mus_df[,-1]
  sigmas<-sigmas_df[,-1]
  cnv_state<-expr
  
  ### hmrf by segmentations
  for (k in c(1:(n_break-1))){
    starts<-breaks_union[k]
    ends<-breaks_union[k+1]-1
    data<-expr[starts:ends,] ### bin * cell
    ###initial
    emission_prob_seg <- emission_prob_mat[starts:ends, , ]
    #print(c(starts,ends))
    R <- cnv_state_init[starts:ends, ]
    #print(paste("starts:", starts, "ends:", ends))
    #print(dim(cnv_state_init))

    #print(R)
    R_last <- R + 1
    
    if (beta_fixed=='au'){
      beta <- get_beta_seg(R, neigh_list)
      if(beta==-1){beta <- beta_default}
      print(beta)
    }else{
      beta <- beta_fixed
    }
    
    print(beta)
    #### estimate the cnv state
    i <- 1
    while(i<=max_iter & sum(R_last[1,]==R[1,])!=n_spot){
      ## E step: update hidden state
      order_idx <- sample(c(1:n_spot), n_spot, F)
      R_last <- R
      ## for each cell in the seg
      for (c in order_idx){
        neigh_ <- neigh_list[[c]]
        x <- data[,c]
        p <- c(1:n_state)
        for (state in 1:n_state){
          potts_prob <- exp(sum(R[,neigh_]==state)*beta)
          p[state] <- prod(emission_prob_seg[,c,state])*potts_prob
        }
        R[,c] <- which.max(p)
      }
      i <- i+1
      
      ## M step: update parameters (beta and gaussian parameters)
      if(beta_fixed=='au'){
        beta_ <- get_beta(R[1,], neigh_list)
        if(beta_!=-1){beta <- beta_}
        print(beta)
      }else{
        beta=beta_fixed
      }
    }
    cnv_state[starts:ends,]<-R
    # print(paste0('inter ',i))
    print(paste0("#### ",k,'/',n_break-1))
  }
  t2 <- Sys.time()
  print(paste0('used ',difftime(t2, t1, units = "mins"),' mins to do hmrf'))
  return(cnv_state)
}


calculate_CNV<- function(norm_count, baseline, 
                         method="BICseq", lambda=0.01,
                         min_interval=2, 
                         smooth=F,
                         sample_dir,
                         maxi_iter=1,
                         lower=60,
                         upper=100,
                         BICseq_dir
){
  bins <- str2bin(colnames(norm_count))
  chr_bkp <- bin2chrbkp(bins)
  if(nrow(norm_count)>1000){
    norm_count<-norm_count[sample(c(1:nrow(norm_count)),1000),]
  }
  
  ## find breakpoints
  print("segmentation")
  BICseq <- BICseq_dir
  lambda <- lambda
  tmp_dir <- sample_dir
  for (q in c(1:maxi_iter)){
    bkp <- c()
    k <- 0
    for(i in unique(bins$chr)){
      print(i)
      if(sum(bins$chr==i)<=2){
        next
      }
      norm_count_ <- norm_count[,bins$chr==i]
      bic_input <- matrix(ncol=2*nrow(norm_count_)+2, 
                          nrow=ncol(norm_count_))
      bic_input[,1:2] <- as.matrix(str2bin(colnames(norm_count_))[,2:3])
      bic_input[,2*c(1:nrow(norm_count_))+1] <- t(norm_count_)
      bic_input[,2*c(1:nrow(norm_count_))+2] <- 
        matrix(data=baseline[bins$chr==i], nrow=ncol(norm_count_), ncol=nrow(norm_count_))
      s <- c(1:ncol(bic_input))
      s[1:2] <- c("start", "end")
      s[2*c(1:nrow(norm_count_))+1] <- 
        rownames(norm_count_)
      s[2*c(1:nrow(norm_count_))+2] <- 
        paste0(rownames(norm_count_),"_ref")
      colnames(bic_input) <- s
      
      write.table(x=bic_input, file=paste0(tmp_dir, i, ".bin"), 
                  quote=F, sep="\t",
                  row.names=F, col.names=T, append=F)
      system(paste0(BICseq, " -i ", paste0(tmp_dir, i, ".bin"),
                    " -l ", lambda))
      seg_re <- read.table(paste0(tmp_dir, i, ".bin_seg"), 
                           header=T, sep="\t")
      bkp_ <- cumsum(seg_re$binNum)
      
      bkp <- merge_bkp(bkp_+k, bkp, min_interval=min_interval)
      k <- k+sum(bins$chr==i)
    }
    for(i in unique(bins$chr)){
      file.remove(paste0(tmp_dir, i, ".bin"))
      file.remove(paste0(tmp_dir, i, ".bin_seg"))
    }
    bkp <- merge_bkp(chr_bkp, bkp, min_interval=2)
    bkp <- sort(bkp)
    #print(paste0("number of bkp: ", length(bkp)))
    if (length(bkp)>=lower &length(bkp)<=upper){
      print(paste0("number of bkp: ", length(bkp)))
      break
    }else if (length(bkp)<lower){
      print("too few breakpoints, changing lambda")
      lambda <- lambda*0.9
    }else{
      print("too many breakpoints, changing lambda")
      lambda <- lambda*1.1
    }
  }
  
  bk_bic<-bkp
  bk_bic[1]<-1
  len_bk<-length(bk_bic)
  if(bk_bic[len_bk]!=ncol(norm_count)){
    bk_bic<-c(bk_bic,ncol(norm_count))
  }
  print(paste0("number of bkp is : ",length(bk_bic)))
  ## save results
  print("saving results")
  saveRDS(bk_bic, file=paste0(sample_dir, "bk_bic.rds"))
  return(bk_bic)
}

state2ratio <- function(cnv_state, copy_ratios){
  tmp <- cnv_state
  for(state in 1:length(copy_ratios)){
    tmp[tmp==state] <- copy_ratios[state]
  }
  return(tmp)
}