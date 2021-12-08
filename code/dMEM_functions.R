############################ Functions used in simulations ############################

library(changepoint)

# function for dMEM with naive kmeans clustering
kmeans_baseline <- function(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,perfect,k=10){
  
  ###############################################################
  # primary: primary source information c(mean,sd,N)
  # supple_gen_mean: vector of generated supplementary means
  # supple_gen_sd: vector of generated supplementary sds
  # supple_gen_N: vector of generated supplementary sample sizes
  # perfect: T as known change point, F as unknown change point
  # k: number of clusters
  ###############################################################
  
  # ground truth of number of exchangeable sources
  num_homo <- length(which(supple_gen_mean == 0))
  supple <- list()
  if(perfect == F){
    # exchangeability unknown
    
    # sample supplementary sources
    for (num_cluster in 1:length(supple_gen_mean)){supple[[num_cluster]] <- c(rnorm(supple_gen_N[num_cluster],supple_gen_mean[num_cluster],supple_gen_sd[num_cluster]))}
    # sample mean
    supple_mean <- unlist(lapply(supple,mean))
    # kmeans cluster supplementary sources
    clustering <- kmeans(supple_mean,k)$cluster
    
  } else {
    # exchangeability known
    
    # sample supplementary sources
    for (num_cluster in 1:num_homo){supple[[num_cluster]] <- c(rnorm(supple_gen_N[num_cluster],supple_gen_mean[num_cluster],supple_gen_sd[num_cluster]))}
    # sample mean
    supple_mean <- unlist(lapply(supple,mean))
    # kmeans cluster supplementary sources
    if(num_homo<=k){
      clustering <- c(1:num_homo)
    } else {
      clustering <- kmeans(supple_mean,k)$cluster
    }
  }
  
  # calculate cluster-level mean/sd/N
  supple_mean_cluster <- c()
  supple_sd_cluster <- c()
  supple_N_cluster <- c()
  for (w in 1:max(clustering)){
    supple_mean_cluster <- c(supple_mean_cluster, mean(unlist(supple[which(clustering==w)])))
    supple_sd_cluster <- c(supple_sd_cluster, sd(unlist(supple[which(clustering==w)])))
    supple_N_cluster <- c(supple_N_cluster, length(unlist(supple[which(clustering==w)])))
  }
  
  # MEM result
  mem_result <- mem_calc(c(mean(primary),sd(primary),length(primary)),
                         supple_mean_cluster,
                         supple_sd_cluster,
                         supple_N_cluster,
                         prior='pi_e')
  postmean <- mem_result$postmean
  postvar <- mem_result$postvar
  postess <- mem_result$ess
  
  return(c('postmean'=postmean,'postvar'=postvar,'ess'=postess))
}

clustering_methods <- function(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,mtd,perfect,badtimes=1,max_num_cluster=10){
  
  ###############################################################
  # primary: primary source information c(mean,sd,N)
  # supple_gen_mean: vector of generated supplementary means
  # supple_gen_sd: vector of generated supplementary sds
  # supple_gen_N: vector of generated supplementary sample sizes
  # mtd: 'bsl' as iMEM and non-random clustering methods ('topk', 'even', 'split')
  # perfect: T as known change point, F as unknown change point
  # badtimes: adding ~ times of extra sources based on marginal weights to the selected supplementary sources
  # max_num_cluster: maximum number of clusters
  ###############################################################
  
  # ground truth of number of exchangeable sources
  num_homo <- length(which(supple_gen_mean == 0))
  
  # generate supplementary sources
  supple <- list()
  for (k in 1:length(supple_gen_mean)){
    supple[[k]] <- c(rnorm(supple_gen_N[k],supple_gen_mean[k],supple_gen_sd[k]))
  }
  
  # sample mean/sd/N
  supple_mean <- unlist(lapply(supple,mean))
  supple_sd <- unlist(lapply(supple,sd))
  supple_N <- supple_gen_N
  
  # marginal MEM
  margmems <- lapply(1:length(supple_N),function(x) {mem_calc(prim=c(mean(primary),sd(primary),length(primary))
                                                              ,means=supple_mean[x],sds=supple_sd[x],
                                                              Ns=supple_N[x],prior='pi_e')} )
  margscores <- sapply(1:length(margmems), function(x){ margmems[[x]]$memlist$postwts[2] })
  
  # changepoint detector
  if (perfect == F){
    # exchangeability unknown
    
    change_point_model <- cpt.mean(sort(margscores,decreasing = T), penalty = "None", method = 'AMOC')
    change_point <- cpts(change_point_model)*badtimes
    
    # source selection
    supple_mean <- supple_mean[order(margscores,decreasing = T)[1:change_point]]
    supple_sd <- supple_sd[order(margscores,decreasing = T)[1:change_point]]
    supple_N <- supple_N[order(margscores,decreasing = T)[1:change_point]]
    supple <- supple[order(margscores,decreasing = T)[1:change_point]]
  } else {
    #exchangeability known
    
    change_point <- num_homo*badtimes
    
    # source selection
    supple_mean <- supple_mean[1:change_point]
    supple_sd <- supple_sd[1:change_point]
    supple_N <- supple_N[1:change_point]
    supple <- supple[1:change_point]
  }
  
  # k: number of sources within each cluster
  # k_split: number of sources within each cluster for the first half of the split case
  # num_cluster: number of clusters
  # num_k_cluster: number of clusters with k sources (the rest of the clusters will have k-1 sources)  
  # num_k_cluster_split: number of clusters with k_split sources for the first half of the split case
  if (change_point < max_num_cluster){
    k <- 1
    num_cluster <- change_point
    num_k_cluster <- change_point
    
    k_split <- 1
    num_k_cluster_split <- max((change_point-round(num_cluster/2)),0)
  } else {
    num_cluster <- max_num_cluster
    k <- ceiling(change_point/num_cluster)
    
    if ((change_point/num_cluster)%%1==0) {
      num_k_cluster <- num_cluster
    } else {
      num_k_cluster <- change_point %% num_cluster
    }
    
    k_split <- round((change_point-round(num_cluster/2)) / (num_cluster-round(num_cluster/2)))
    if (((change_point-round(num_cluster/2)) / (num_cluster-round(num_cluster/2)))%%1==0) {
      num_k_cluster_split <- num_cluster-round(num_cluster/2)
    } else {
      num_k_cluster_split <- (change_point-round(num_cluster/2)) %% (num_cluster-round(num_cluster/2))
    }
    }
    
  supple_mean_cluster <- c()
  supple_sd_cluster <- c()
  supple_N_cluster <- c()
  
  if(mtd=='bsl'){ 
    # iMEM
    
    mem_result <- imem_marg(c(mean(primary),sd(primary),length(primary)),
                            supple_mean,supple_sd,supple_N,
                            prior='pi_e',final_grpsize = min(max_num_cluster,num_cluster))
    
  } else if (mtd=='topk') {
    # topk clustering
    
    # calculate cluster-level mean/sd/N
    for (j in 1:num_cluster){
      if (j <= num_k_cluster) {
        supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[((j-1)*k+1):min(j*k,length(supple))])))
        supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[((j-1)*k+1):min(j*k,length(supple))])))
        supple_N_cluster <- c(supple_N_cluster, sum(supple_N[((j-1)*k+1):min(j*k,length(supple))]))
      } else {
        supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                                          min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))])))
        supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                                    min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))])))
        supple_N_cluster <- c(supple_N_cluster, sum(supple_N[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                               min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))]))
      }
    }
    
    supple_mean_cluster <- supple_mean_cluster[is.na(supple_N_cluster)==F]
    supple_sd_cluster <- supple_sd_cluster[is.na(supple_N_cluster)==F]
    supple_N_cluster <- supple_N_cluster[is.na(supple_N_cluster)==F]
    
    # MEM result
    mem_result <- mem_calc(c(mean(primary),sd(primary),length(primary)),
                           supple_mean_cluster,
                           supple_sd_cluster,
                           supple_N_cluster,
                           prior='pi_e')
    
  } else if (mtd=='even') {
    # evenly distributed clustering
    
    # calculate cluster-level mean/sd/N
    for (j in 1:num_cluster){
      index <- seq(1,length(supple),num_cluster)+j-1
      supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[index[which(index <= length(supple))]])))
      supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[index[which(index <= length(supple))]])))
      supple_N_cluster <- c(supple_N_cluster, sum(supple_N[index[which(index <= length(supple))]]))
    }
    
    # MEM result
    mem_result <- mem_calc(c(mean(primary),sd(primary),length(primary)),
                           supple_mean_cluster,
                           supple_sd_cluster,
                           supple_N_cluster,
                           prior='pi_e')
    
  } else if (mtd=='split') {
    # top half & combine half clustering
    
    # calculate cluster-level mean/sd/N
    for (j in 1:round(num_cluster/2)){
      supple_mean_cluster <- c(supple_mean_cluster,supple_mean[j])
      supple_sd_cluster <- c(supple_sd_cluster,supple_sd[j])
      supple_N_cluster <- c(supple_N_cluster, supple_N[j])
    }
    for (j in 1:(num_cluster-round(num_cluster/2))){
      if (j <= num_k_cluster_split){
        supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[((j-1)*k_split+round(num_cluster/2)+1):
                                                                          min(j*k_split+round(num_cluster/2),length(supple))])))
        supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[((j-1)*k_split+round(num_cluster/2)+1):
                                                                    min(j*k_split+round(num_cluster/2),length(supple))])))
        supple_N_cluster <- c(supple_N_cluster, sum(supple_N[((j-1)*k_split+round(num_cluster/2)+1):
                                                               min(j*k_split+round(num_cluster/2),length(supple))]))
      } else {
        supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[(num_k_cluster_split*k+(j-num_k_cluster_split-1)*(k-1)+round(num_cluster/2)+1):
                                                                          min(num_k_cluster_split*k+(j-num_k_cluster_split)*(k-1)+round(num_cluster/2),length(supple))])))
        supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[(num_k_cluster_split*k+(j-num_k_cluster_split-1)*(k-1)+round(num_cluster/2)+1):
                                                                    min(num_k_cluster_split*k+(j-num_k_cluster_split)*(k-1)+round(num_cluster/2),length(supple))])))
        supple_N_cluster <- c(supple_N_cluster, sum(supple_N[(num_k_cluster_split*k+(j-num_k_cluster_split-1)*(k-1)+round(num_cluster/2)+1):
                                                               min(num_k_cluster_split*k+(j-num_k_cluster_split)*(k-1)+round(num_cluster/2),length(supple))]))
      }

    }
    
    # MEM result
    mem_result <- mem_calc(c(mean(primary),sd(primary),length(primary)),
                           supple_mean_cluster,
                           supple_sd_cluster,
                           supple_N_cluster,
                           prior='pi_e')
    
  }
  
  postmean <- mem_result$postmean
  postvar <- mem_result$postvar
  postess <- mem_result$ess
  
  return(c('postmean'=postmean,'postvar'=postvar,'ess'=postess))
}

random_clustering <- function(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,m,perfect,badtimes=1,max_num_cluster=10){
  
  ###############################################################
  # primary: primary source information c(mean,sd,N)
  # supple_gen_mean: vector of generated supplementary means
  # supple_gen_sd: vector of generated supplementary sds
  # supple_gen_N: vector of generated supplementary sample sizes
  # m: run m times of random clustering and take the average of the results
  # perfect: T as known change point, F as unknown change point
  # badtimes: adding ~ times of extra sources based on marginal weights to the selected supplementary sources
  # max_num_cluster: maximum number of clusters
  ###############################################################
  
  # ground truth of number of exchangeable sources
  num_homo <- length(which(supple_gen_mean == 0))
  
  # generate supplementary sources
  supple <- list()
  for (k in 1:length(supple_gen_mean)){
    supple[[k]] <- c(rnorm(supple_gen_N[k],supple_gen_mean[k],supple_gen_sd[k]))
  }
  
  # sample mean/sd/N
  supple_mean <- unlist(lapply(supple,mean))
  supple_sd <- unlist(lapply(supple,sd))
  supple_N <- supple_gen_N
  
  # marginal MEM
  margmems <- lapply(1:length(supple_N),function(x) {mem_calc(prim=c(mean(primary),sd(primary),length(primary))
                                                              ,means=supple_mean[x],sds=supple_sd[x],
                                                              Ns=supple_N[x],prior='pi_e')} )
  margscores <- sapply(1:length(margmems), function(x){ margmems[[x]]$memlist$postwts[2] })
  
  # changepoint detector
  if (perfect == F){
    # changepoint unknown
    
    change_point_model <- cpt.mean(sort(margscores,decreasing = T), penalty = "None", method = 'AMOC')
    change_point <- cpts(change_point_model)*badtimes
    
    supple_mean <- supple_mean[order(margscores,decreasing = T)[1:change_point]]
    supple_sd <- supple_sd[order(margscores,decreasing = T)[1:change_point]]
    supple_N <- supple_N[order(margscores,decreasing = T)[1:change_point]]
    supple <- supple[order(margscores,decreasing = T)[1:change_point]]
  } else {
    # changepoint known
    
    change_point <- num_homo*badtimes
    
    supple_mean <- supple_mean[1:change_point]
    supple_sd <- supple_sd[1:change_point]
    supple_N <- supple_N[1:change_point]
    supple <- supple[1:change_point]
  }
  
  # k: number of sources within each cluster
  # num_cluster: number of clusters
  # num_k_cluster: number of clusters with k sources (the rest of the clusters will have k-1 sources)
  if (change_point < max_num_cluster){
    k <- 1
    num_cluster <- change_point
    num_k_cluster <- change_point

  } else {
    num_cluster <- max_num_cluster
    k <- ceiling(change_point/num_cluster)
    if ((change_point/num_cluster)%%1==0) {
      num_k_cluster <- num_cluster
    } else {
      num_k_cluster <- change_point %% num_cluster
    }

  }
  
  cluster_mem_postmean <- c()
  cluster_mem_postvar <- c()
  cluster_mem_ess <- c()
  
  for (b in 1:m) {
    # random clustering m times
    
    shuffle <- sample(c(1:length(supple)))
    # shuffle the selected groups
    supple_mean <- supple_mean[shuffle]
    supple_sd <- supple_sd[shuffle]
    supple_N <- supple_N[shuffle]
    supple <- supple[shuffle]
    
    # calculate cluster-level mean/sd/N
    supple_mean_cluster <- c()
    supple_sd_cluster <- c()
    supple_N_cluster <- c()
    for (j in 1:num_cluster){
      if (j <= num_k_cluster) {
        supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[((j-1)*k+1):min(j*k,length(supple))])))
        supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[((j-1)*k+1):min(j*k,length(supple))])))
        supple_N_cluster <- c(supple_N_cluster, sum(supple_N[((j-1)*k+1):min(j*k,length(supple))]))
      } else {
        supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                                          min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))])))
        supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                                    min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))])))
        supple_N_cluster <- c(supple_N_cluster, sum(supple_N[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                               min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))]))
      }
    }
    
    supple_mean_cluster <- supple_mean_cluster[is.na(supple_N_cluster)==F]
    supple_sd_cluster <- supple_sd_cluster[is.na(supple_N_cluster)==F]
    supple_N_cluster <- supple_N_cluster[is.na(supple_N_cluster)==F]
    
    # MEM result
    mem_result <- mem_calc(c(mean(primary),sd(primary),length(primary)),
                           supple_mean_cluster,
                           supple_sd_cluster,
                           supple_N_cluster,
                           prior='pi_e')
    
    # record results of the random clustering
    cluster_mem_postmean <- c(cluster_mem_postmean, mem_result$postmean)   
    cluster_mem_postvar <- c(cluster_mem_postvar, mem_result$postvar)   
    cluster_mem_ess <- c(cluster_mem_ess, mem_result$ess)   
  }
  
  # simple averaging
  postmean <- mean(cluster_mem_postmean)
  postvar <- mean(cluster_mem_postvar)
  postess <- mean(cluster_mem_ess)
  
  return(c('postmean'=postmean,'postvar'=postvar,'ess'=postess))
}


############################ dMEM functions used in application ############################

credint <- function(memobj, cred.lev=0.95){
  nsamp.mod <- c(rmultinom(1,1000000,memobj$memlist$postwts))
  
  post.dist <- unlist(apply(cbind(nsamp.mod,memobj$memlist$pmeans,memobj$memlist$pvars), 1, function(x){rnorm(x[1],x[2],x[3])}))
  
  hpdint <- boa.hpd(post.dist,alpha=1-cred.lev)
  upr <- hpdint[1]
  lwr <- hpdint[2]
  
  return(c(upr,lwr))
}

DMEM_app <- function(primary,supple_mean,supple_sd,supple_N,supple,mtd='random',m=1,max_num_cluster=10){
  
  ###############################################################
  # primary: primary source information c(mean,sd,N)
  # supple_gen_mean: vector of generated supplementary means
  # supple_gen_sd: vector of generated supplementary sds
  # supple_gen_N: vector of generated supplementary sample sizes
  # supple: sampled supplementary sources
  # mtd: clustering methods ('topk', 'even', 'random')
  # m: run m times of random clustering and take the average of the results
  # max_num_cluster: maximum number of clusters
  ###############################################################
  
  # marginal MEM
  margmems <- lapply(1:length(supple_N),function(x) {mem_calc(prim=primary
                                                              ,means=supple_mean[x],sds=supple_sd[x],
                                                              Ns=supple_N[x],prior='pi_e')} )
  margscores <- sapply(1:length(margmems), function(x){ margmems[[x]]$memlist$postwts[2] })
  
  # detect change point
  change_point_model <- cpt.mean(sort(margscores,decreasing = T), penalty = "None", method = 'AMOC')
  change_point <- cpts(change_point_model)
  
  # hard threshold
  hard_threshold <- 0
  cutoff <- max(which(sort(margscores,decreasing = T) >= hard_threshold),1)
  if (cutoff < change_point){
    change_point <- cutoff
  }
  
  # source selection
  supple_mean <- supple_mean[order(margscores,decreasing = T)[1:change_point]]
  supple_sd <- supple_sd[order(margscores,decreasing = T)[1:change_point]]
  supple_N <- supple_N[order(margscores,decreasing = T)[1:change_point]]
  supple <- supple[order(margscores,decreasing = T)[1:change_point]]
  
  # k: number of sources within each cluster
  # num_cluster: number of clusters
  # num_k_cluster: number of clusters with k sources (the rest of the clusters will have k-1 sources)
  if (change_point < max_num_cluster){
    k <- 1
    num_cluster <- change_point
    num_k_cluster <- change_point
    
  } else {
    num_cluster <- max_num_cluster
    k <- ceiling(change_point/num_cluster)
    if ((change_point/num_cluster)%%1==0) {
      num_k_cluster <- num_cluster
    } else {
      num_k_cluster <- change_point %% num_cluster
    }
    
  }
  
  # calculate cluster-level mean/sd/N
  supple_mean_cluster <- c()
  supple_sd_cluster <- c()
  supple_N_cluster <- c()
  
  if (mtd=='random'){  
    # random clustering or random averaging
    
    cluster_mem_postmean <- c()
    cluster_mem_postvar <- c()
    cluster_mem_ess <- c()
    cluster_mem_upr <- c()
    cluster_mem_lwr <- c()
    for (b in 1:m) {
      
      # shuffle the selected sources
      shuffle <- sample(c(1:length(supple)))
      supple_mean <- supple_mean[shuffle]
      supple_sd <- supple_sd[shuffle]
      supple_N <- supple_N[shuffle]
      supple <- supple[shuffle]
      
      supple_mean_cluster <- c()
      supple_sd_cluster <- c()
      supple_N_cluster <- c()
      
      for (j in 1:num_cluster){
        if (j <= num_k_cluster) {
          supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[((j-1)*k+1):min(j*k,length(supple))])))
          supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[((j-1)*k+1):min(j*k,length(supple))])))
          supple_N_cluster <- c(supple_N_cluster, sum(supple_N[((j-1)*k+1):min(j*k,length(supple))]))
        } else {
          supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                                            min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))])))
          supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                                      min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))])))
          supple_N_cluster <- c(supple_N_cluster, sum(supple_N[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                                 min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))]))
        }
      }
      
      supple_mean_cluster <- supple_mean_cluster[is.na(supple_N_cluster)==F]
      supple_sd_cluster <- supple_sd_cluster[is.na(supple_N_cluster)==F]
      supple_N_cluster <- supple_N_cluster[is.na(supple_N_cluster)==F]
      
      # MEM result
      mem_result <- mem_calc(primary,
                             supple_mean_cluster,
                             supple_sd_cluster,
                             supple_N_cluster,
                             prior='pi_e')
      
      cluster_mem_postmean <- c(cluster_mem_postmean, mem_result$postmean)   
      cluster_mem_postvar <- c(cluster_mem_postvar, mem_result$postvar)   
      cluster_mem_ess <- c(cluster_mem_ess, mem_result$ess)   
      interval <- credint(mem_result)
      cluster_mem_upr <- c(cluster_mem_upr, interval[1])
      cluster_mem_lwr <- c(cluster_mem_lwr, interval[2])
      
    }
    
    # simple averaging
    postmean <- mean(cluster_mem_postmean)
    postvar <- mean(cluster_mem_postvar)
    postess <- mean(cluster_mem_ess)  
    postupr <- mean(cluster_mem_upr)
    postlwr <- mean(cluster_mem_lwr)
  
  } else if (mtd=='topk') {
    # topk clustering
    
    for (j in 1:num_cluster){
      if (j <= num_k_cluster) {
        supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[((j-1)*k+1):min(j*k,length(supple))])))
        supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[((j-1)*k+1):min(j*k,length(supple))])))
        supple_N_cluster <- c(supple_N_cluster, sum(supple_N[((j-1)*k+1):min(j*k,length(supple))]))
      } else {
        supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                                          min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))])))
        supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                                    min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))])))
        supple_N_cluster <- c(supple_N_cluster, sum(supple_N[(num_k_cluster*k+(j-num_k_cluster-1)*(k-1)+1):
                                                               min(num_k_cluster*k+(j-num_k_cluster)*(k-1),length(supple))]))
      }
    }
    
    supple_mean_cluster <- supple_mean_cluster[is.na(supple_N_cluster)==F]
    supple_sd_cluster <- supple_sd_cluster[is.na(supple_N_cluster)==F]
    supple_N_cluster <- supple_N_cluster[is.na(supple_N_cluster)==F]
    
    # MEM result
    mem_result <- mem_calc(primary,
                           supple_mean_cluster,
                           supple_sd_cluster,
                           supple_N_cluster,
                           prior='pi_e')
    postmean <- mem_result$postmean
    postvar <- mem_result$postvar
    postess <- mem_result$ess
    interval <- credint(mem_result)
    postupr <- interval[1]
    postlwr <- interval[2]
    
  } else if (mtd=='even') {
    # evenly distributed clustering
    
    for (j in 1:num_cluster){
      index <- seq(1,length(supple),num_cluster)+j-1
      supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[index[which(index <= length(supple))]])))
      supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[index[which(index <= length(supple))]])))
      supple_N_cluster <- c(supple_N_cluster, sum(supple_N[index[which(index <= length(supple))]]))
    }
    
    # MEM result
    mem_result <- mem_calc(primary,
                           supple_mean_cluster,
                           supple_sd_cluster,
                           supple_N_cluster,
                           prior='pi_e')
    postmean <- mem_result$postmean
    postvar <- mem_result$postvar
    postess <- mem_result$ess
    interval <- credint(mem_result)
    postupr <- interval[1]
    postlwr <- interval[2]
    
  }
  
  return(c('postmean'=postmean,'postvar'=postvar,'ess'=postess,'upr'=postupr,'lwr'=postlwr))
}

