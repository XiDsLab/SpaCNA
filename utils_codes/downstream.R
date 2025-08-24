plot_heatmap <- function(copy_ratio,
                         del_chrXY=T,
                         bin="none",
                         cell_cluster="none",
                         cluster_method="kmeans", K=5,
                         cell_annotation="none",
                         save_hc=F,
                         output_dir="./",
                         output_name="copy_ratio.png"
){
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  save_dir <- paste0(output_dir, output_name)
  if(!identical(bin, "none")){
    bins <- bin[,2:5]
  }else{
    bins <- str2bin(colnames(copy_ratio))
  }
  
  if(del_chrXY){
    f <- !(bins$chr %in% c("chrX","chrY","X","Y"))
    copy_ratio <- copy_ratio[,f]
    bins <- bins[f,]
  }
  chr_bkp <- bin2chrbkp(bins)
  
  ## cluster
  print("clustering")
  if(min(copy_ratio)<=0){
    temp <- copy_ratio
  }else{
    temp <- log(copy_ratio+0.01)
  }
  temp <- copy_ratio
  
  if(identical(cell_cluster,"none")){
    d <- parallelDist::parallelDist(temp)
    hc_re <- hclust(d, method="ward.D2")
  }else if(identical(cell_cluster,"kmeans")){
    kmeans_pc <- 50
    if(identical(kmeans_pc, "none")){
      temp_ <- temp
    }else{
      temp_ <- irlba::prcomp_irlba(temp, n=kmeans_pc)$x
    }
    group <- kmeans(temp_, centers=K)$cluster
    hc_re <- ComplexHeatmap::cluster_within_group(t(temp), group)
  }else{
    cell_cluster <- as.numeric(as.factor(cell_cluster))
    hc_re <- ComplexHeatmap::cluster_within_group(t(temp), cell_cluster)
    K <- length(unique(cell_cluster))
  }
  
  ## plot
  print("plotting")
  plot_data <- copy_ratio
  
  color_thre <- c(0.5, 1.5)
  c <- colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")))(100)
  colormap <-
    circlize::colorRamp2(
      c(1:100)*(color_thre[2]-color_thre[1])*0.01+color_thre[1], c)
  
  ## set row annotation
  if(identical(cell_annotation,"none")){
    a <- ComplexHeatmap::HeatmapAnnotation(
      cell_type=ComplexHeatmap::anno_empty(border=F),which="row")
  }else if(all(cell_annotation %in% c("tumor","normal"))){
    a <- ComplexHeatmap::HeatmapAnnotation(
      cell_type=cell_annotation,
      col=list(cell_type=c("normal"="gray",
                           "tumor"="purple")),
      which="row",show_annotation_name=F)
  }else{
    a <- ComplexHeatmap::HeatmapAnnotation(
      cell_type=cell_annotation, which="row",show_annotation_name=F)
  }
  
  ## set col annotation
  n_chr <- length(unique(bins$chr))
  temp <- rep("gray",n_chr)
  temp[2*c(1:floor(n_chr/2))] <- "white"
  names(temp) <- unique(bins$chr)
  b <- ComplexHeatmap::HeatmapAnnotation(
    chr=bins$chr, col=list(chr=temp), which="column", show_annotation_name=F)
  
  ## heatmap
  # ComplexHeatmap::ht_opt$message <- FALSE
  ht <- ComplexHeatmap::Heatmap(
    plot_data, name="copy ratio", col=colormap,
    cluster_rows=hc_re, cluster_columns=FALSE,
    show_row_names=FALSE, show_column_names=FALSE,
    split=K,
    left_annotation=a, top_annotation=b)
  
  ## add vertical lines
  k <- K
  if(substr(save_dir, nchar(save_dir)-2, nchar(save_dir))=="png"){
    png(filename=save_dir, width=2000, height=1000)
  }else{
    pdf(file=save_dir, width=20, height=10)
  }
  print(ht)
  for(i in 1:k){
    for(j in chr_bkp){
      ComplexHeatmap::decorate_heatmap_body(
        "copy ratio",
        {grid::grid.lines(c(j/dim(plot_data)[2], j/dim(plot_data)[2]),
                          c(0,1), gp=grid::gpar(lty=1, lwd=2))},
        slice=i)
    }
  }
  dev.off()
  
  ## save cluster result
  if(save_hc){
    saveRDS(hc_re, paste0(output_dir, "heatmap_hc.rds"))
  }
}
