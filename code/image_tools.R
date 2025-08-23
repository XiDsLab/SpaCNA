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

get_spot_neighbors <- function(spot_location, k=7, dist_thre=500){
  neigh_list<-list()
  knn<- FindNeighbors(as.matrix(spot_location), 
                      k.param=k, 
                      annoy.metric="manhattan", 
                      return.neighbor=TRUE)
  for(i in 1:nrow(spot_location)){
    idx1 <- knn@nn.idx[i,-1] 
    dist1 <- knn@nn.dist[i,-1]
    idx1 <- idx1[dist1 <= dist_thre] 
    neigh_list[[i]]<- idx1
  }
  names(neigh_list) <- rownames(spot_location)
  return(neigh_list)
}


