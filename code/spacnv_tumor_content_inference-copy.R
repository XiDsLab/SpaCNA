source("./spatial_cnv/code/code_ver2/utils.R")
source("./spatial_cnv/code/code_ver2/downstream.R")
source("./spatial_cnv/code/code_ver2/normalize.R")
source("./spatial_cnv/code/code_ver2/hrmf_em.R")
source("./spatial_cnv/code/code_ver2/hrmf_init.R")

sample_dir <- "./spatial_cnv/data/breast/"
plot_dir <- "./spatial_cnv/bc_analysis/brca1/"
obj <- readRDS(paste0(sample_dir, "seurat_object.rds"))
## input
cns <- readRDS(paste0(sample_dir, "cns.rds")) # copy number state
count_norm <- readRDS(paste0(sample_dir, "count_norm.rds")) # copy ratio
bk_bic <- readRDS(paste0(sample_dir, "bk_bic.rds")) # breakpoints
cell_barcodes <- colnames(count_norm)
K <- 7 # number of cluster
cnv_thre <- 8 # number of minimum cnv segment 

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

a<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/tumor4_perc.rds")
b<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/tumor4_perc_esti2.rds")
c<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/hrmf_para.rds")

pearson_corr <- cor(a$tumor_perc, b$p_estimate, method = "pearson")
pearson_corr <- cor(a$tumor_perc[181:400], b$p_estimate[181:400], method = "pearson")


# 创建示例数据
data <- data.frame(
  length = c(1:nrow(a), 1:nrow(a)),  # 数据点的索引
  value = c(a$tumor_perc, b$p_estimate),  # x 和 y 的值
  variable = rep(c("true", "estimate"), each = nrow(a))  # 标识 x 和 y 的变量
)

# 创建散点图
ggplot(data, aes(x = length, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("true" = "red", "estimate" = "green")) +
  ggtitle("Scatter Plot of x and y") +
  xlab("Length") +
  ylab("Value") +
  theme_minimal()

a2<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/copy_ratio_purei.rds")
a3<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/copy_ratio_segi.rds")
f<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/f.rds")
cns<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/copy_ratio_sim_hmrf_bic_ver2.rds")
count_norm<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/count_norm.rds")
count_base<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/count_base_filter.rds")
idx<-rownames(count_norm)[f]
sel<-count_base[idx,"tumor2_22"]
a2[which(sel==1)]
a3[which(sel==1)]
df_2<-data.frame(x=a2[which(sel!=1)]- 1,y=a3[which(sel!=1)]-1)
plot(df_2$x, df_2$y, xlab = "Copy Ratio Pure", ylab = "Copy Ratio Seg", pch = 19)

idx2<-(a2-1<0&a3-1>0)
idx2<-idx[idx2]
count_base[idx2,"tumor2_22"]
sel2<-count_base[idx2,"tumor2_22"]
count_norm[idx2,"tumor2_22"]
data <- data.frame(x = a$tumor_per, y = b$p_estimate)
sel2<-cns[idx,"tumor2_22"]
df2<-data.frame(true1=sel,esti=sel2)
# 创建相关性散点图
ggplot(data, aes(x = x, y = y)) +
  geom_point(color = "blue", size = 3) +  # 绘制散点
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # 添加线性回归线
  ggtitle("Scatter Plot of x and y with Regression Line") +  # 设置标题
  xlab("x values") +  # 设置 x 轴标签
  ylab("y values") +  # 设置 y 轴标签
  theme_minimal()  # 使用最小主题

a<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/tumor2_perc.rds")
a1<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/copy_ratio_pure_gain_names.rds")
a2<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/copy_ratio_pure_loss.rds")
sort(a$tumor_perc,decreasing = T)[1:10]
idx<-(a$tumor_perc %in% sort(a$tumor_perc,decreasing = T)[1:10])
rownames(a)[idx]
pearson_corr <- cor(a$tumor_perc[181:400], b$p_estimate[181:400], method = "pearson")
hmrf<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/copy_ratio_sim_hmrf_bic_ver2.rds")
count_norm<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/count_norm2.rds")
test<-count_norm["chr6_3.7e+07_3.8e+07",181:400]
plot(test)
test2<-hmrf["chr6_3.7e+07_3.8e+07",181:400]
tumor_content<-a$tumor_perc[181:400]
names(tumor_content)<-rownames(a)[181:400]
aa<-sort(tumor_content)
test<-test[names(aa)]
test2<-test2[names(aa)]
plot(test)

df_<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/df_.rds")
#count_norm<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/count_norm.rds")
#test1<-count_norm[,]
# 假设df_已经创建，并且re是拟合的线性模型
# 画出散点图
plot(df_$x, df_$y, xlab = "Copy Ratio Pure", ylab = "Copy Ratio Seg", pch = 19)

# 在散点图上添加拟合线
abline(re, col = "blue")

# 计算R²值
summary_re <- summary(re)
r_squared <- summary_re$r.squared

# 打印R²值
print(paste("R-squared:", r_squared))

copy_ratio_pure_gain <- matrix(1,nrow=K,ncol= n_seg)
copy_ratio_pure_loss <- matrix(1,nrow=K,ncol= n_seg)
for(k in 1:K){
  clone_cell <- cell_barcodes[clone == k]
  for(i in 1:n_seg){
    tmp <- copy_ratio_seg[clone_cell, i]
    f <- cns_seg[clone_cell,i] < 1
    if(sum(f) > 20){
      copy_ratio_pure_loss[k,i] <- mean(sort(tmp[f])[1:20])
    }
    f <- cns_seg[clone_cell,i] > 1
    if(sum(f) >20){
      copy_ratio_pure_gain[k,i] <- mean(sort(tmp, decreasing = T)[1:20])
    }
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
  copy_ratio_pure[f] <- copy_ratio_pure_gain[k,f]
  f <- clone_cnv < 1
  copy_ratio_pure[f] <- copy_ratio_pure_loss[k,f]
  
  clone_cell <- cell_barcodes[clone == k]
  tumor_content_estimate <- rep(0, length(clone_cell))
  for(i in 1:length(clone_cell)){
    f <- clone_cnv!= 1
    df_ <- data.frame(x = copy_ratio_pure[f] - 1, 
                      y = copy_ratio_seg[clone_cell[i],f] - 1)
    re <- lm(data = df_, formula = y ~ x + 0)
    df[clone_cell[i], "p_estimate"] <- re$coefficients
    df[clone_cell[i], "p_estimate_confidence"] <- summary(re)$adj.r.squared
    if(clone_cell[i]=="tumor2_22"){
      saveRDS(df_,paste0(plot_dir_sim, "df_.rds"))
      saveRDS(re,paste0(plot_dir_sim, "re.rds"))
      saveRDS(f,paste0(plot_dir_sim, "f.rds"))
      saveRDS(copy_ratio_pure[f],paste0(plot_dir_sim, "copy_ratio_purei.rds"))
    }
  }
}

get_corr_aucbybin_accuracy<-function(cnv_mat,true_cnv_mat,tumor_idx){
  
  inter1<-intersect(rownames(cnv_mat), rownames(true_cnv_mat))
  mat1<-cnv_mat[inter1,tumor_idx]
  mat2<-true_cnv_mat[inter1,tumor_idx]
  correlations=c()
  for (i in 1:ncol(mat1)) {
    col1 <- mat1[, i]
    col2 <- mat2[, i]
    correlations[i] <- cor(col1, col2)
  }
  true_cnv<-abs(mat2-1)
  true_cnv[true_cnv!=0]<-1### 0 for no cnv, 1 for cnv
  pre_cnv<-abs(mat1-1)
  auc_bybin <- auc(roc(c(true_cnv), c(pre_cnv)))
  accuracy<-sum(mat1==mat2)/(nrow(mat1)*ncol(mat1))
  return(list(correlations=correlations,auc_bybin=auc_bybin,accuracy=accuracy))
}


get_false_negtive3<-function(cnv_mat,true_cnv_mat,tumor_idx){
  inter1<-intersect(rownames(cnv_mat), rownames(true_cnv_mat))
  mat1<-cnv_mat[inter1,tumor_idx]
  mat2<-true_cnv_mat[inter1,tumor_idx]
  predicted_labels<-c(mat1)
  actual_labels<-c(mat2)
  ###true negative
  neg_idx<-actual_labels==1
  t_n<-sum(predicted_labels[neg_idx]==1)
  
  ###false negativa
  pos_idx<-actual_labels!=1
  f_n<-sum(predicted_labels[pos_idx]==1)
  
  ###true positive
  t_p<-sum(predicted_labels[pos_idx]==actual_labels[pos_idx])
  
  ###false positive
  idx_05<-actual_labels==0.5
  idx_15<-actual_labels==1.5
  print(sum(predicted_labels[idx_05]==1.5)+sum(predicted_labels[idx_15]==0.5))
  print(sum(pos_idx))
  print(sum(predicted_labels[neg_idx]!=1))
  print(sum(neg_idx))
  f_p<-(sum(predicted_labels[idx_05]==1.5)+
          sum(predicted_labels[idx_15]==0.5)+
          sum(predicted_labels[neg_idx]!=1))
  
  #####ratio
  precision<-t_p/(t_p+f_p)
  recall<-t_p/(t_p+f_n)
  f1<-2*precision*recall/(precision+recall)
  
  
  return(
    list(true_negative=t_n,
         false_negative=f_n,
         true_positive=t_p,
         false_positive=f_p,
         precision=precision,
         recall=recall,
         f1=f1)
  )
}

copy_ratio_sim2<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/simulate/copy_ratio_sim_hmrf_bic_ver2.rds")
count_base_filter<-readRDS("/Users/paixiaoshan/Documents/Rdocument/linux/simulate/count_base_filter.rds")
tumor_idx<-colnames(copy_ratio_sim2)
matric_spa2<-get_false_negtive3(copy_ratio_sim2,count_base_filter,tumor_idx)








