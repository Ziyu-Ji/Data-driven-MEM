############################Simulation Scenaio 1 ############################
################# Nonexchangeable sources with the same mean ################
#############################################################################

source('iMEM_functions.R')
source('dMEM_functions.R')

# primary source
true_data <- rnorm(100000,0,sd=1) # The true distribution N(0,sd=1) for primary data source

# supplementary source
sup_sd <- runif(100000,min=0.5,max=1.5) # The distribution of sd for supplementary data sources
sup_N <- c(15:25) # The sample size range for supplementary data sources
sup_homo <- 0 # The homogeneous (exchangeable) supplementry mean
sup_hetero_mean <- 1 # The heterogeneous (nonexchangeable) supplementary mean

# function running the simulations for kmeans, iMEM and dMEM with non-random clustering methods
simulation_homo_cluster <- function(mtd='bsl',perfect=F,k=10,q=10){ 
  
###################################
# mtd: method, including naive kmeans ('kmeans'), iMEM ('bsl'), and non-random clustering methods ('topk', 'even', 'split')
# perfect: use changepoint detector or not
# k: number of clusters for kmeans
# q: number of selected sources for iMEM
###################################
  
  # results are generated from 1000 random simulations
  for (i in 1:1000){
    
    # sampling primary source
    primary <- true_data[sample(1:length(true_data), 20)] # The primary source is sampled from the true data with sample size = 20
    
    # sampling supplementary sources
    supple_gen_mean <- c(rep(sup_homo,60),rep(sup_hetero_mean,40)) # 60 homogeneous & 40 heterogeneous
    supple_gen_sd <- sup_sd[sample(1:100000,100)] 
    supple_gen_N <- sup_N[sample(1:length(sup_N),100,replace = T)]
    
    # results for kmeans, bsl (iMEM), or dMEM
    if (mtd=='kmeans'){
      mem_result <- kmeans_baseline(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,perfect=perfect,k=k)
    } else if(mtd=='bsl'){
      mem_result <- clustering_methods(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,mtd=mtd,perfect=perfect,max_num_cluster=q)
    } else {
      mem_result <- clustering_methods(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,mtd=mtd,perfect=perfect)
    }
    
    # generate result table
    if (i==1){
      df <- t(c(mem_result['postmean'],mem_result['postvar'],mem_result['ess']))
      df <- as.data.frame(df)
      colnames(df) <- c('postmean','postvar','postess')
      
      write.csv(df,paste0('optimal_clustering_simulation_',sup_hetero_mean,'heteromean_',mtd,'method_',perfect,'perfect_',min(k,q),'kq.csv'),row.names=F)
      
    } else {
      df <- t(c(mem_result['postmean'],mem_result['postvar'],mem_result['ess']))
      write.table(df,paste0('optimal_clustering_simulation_',sup_hetero_mean,'heteromean_',mtd,'method_',perfect,'perfect_',min(k,q),'kq.csv'),append=T,sep=',',row.names=F,col.names=F)
    }
  }
  
}


# function running the simulations for dMEM with random clustering methods
simulation_random_cluster <- function(m, perfect=F, max_num_cluster=10){ 
  
  ###################################
  # m: times of randomly clustering within one iteration, results are the average of the m times
  # perfect: use changepoint detector or not
  # max_num_cluster: number of clusters
  ###################################
  
  # results are generated from 1000 random simulations
  for (i in 1:1000){
    
    # sampling primary source
    primary <- true_data[sample(1:10000, 20)] # The primary source is sampled from the true data with sample size = 20
    
    # sampling supplementary sources
    supple_gen_mean <- c(rep(sup_homo,60),rep(sup_hetero_mean,40))
    supple_gen_sd <- sup_sd[sample(1:10000,100)]
    supple_gen_N <- sup_N[sample(1:length(sup_N),100,replace = T)]
    
    # dMEM results for random clustering
    mem_result <- random_clustering(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,m,perfect=perfect,max_num_cluster=max_num_cluster)
    
    # generate result table
    if (i==1){
      df <- t(c(mem_result['postmean'],mem_result['postvar'],mem_result['ess']))
      df <- as.data.frame(df)
      colnames(df) <- c('postmean','postvar','postess')
      
      write.csv(df,paste0('optimal_clustering_simulation_',m,'iter_',sup_hetero_mean,'heteromean_random_',perfect,'perfect_',max_num_cluster,'max_clst.csv'),row.names=F)
      
    } else {
      df <- t(c(mem_result['postmean'],mem_result['postvar'],mem_result['ess']))
      write.table(df,paste0('optimal_clustering_simulation_',m,'iter_',sup_hetero_mean,'heteromean_random_',perfect,'perfect_',max_num_cluster,'max_clst.csv'),append=T,sep=',',row.names=F,col.names=F)
    }
  }
  
}

simulation_homo_cluster(mtd='kmeans',k=5)
simulation_homo_cluster(mtd='kmeans')
simulation_homo_cluster(mtd='bsl',q=5)
simulation_homo_cluster(mtd='bsl')
simulation_homo_cluster(mtd='topk')
simulation_homo_cluster(mtd='even')
simulation_homo_cluster(mtd='split')

simulation_random_cluster(10) # random averaging
simulation_random_cluster(1,perfect=F, max_num_cluster=1)
simulation_random_cluster(1,perfect=F, max_num_cluster=3)
simulation_random_cluster(1,perfect=F, max_num_cluster=10)

############################Simulation Scenaio 2 ############################
################### Varying degrees of nonexchangeability  ##################
#############################################################################

# primary source
true_data <- rnorm(100000,0,sd=1) # The true distribution N(0,sd=1) for primary data source

# supplementary source
sup_sd <- runif(100000,min=0.5,max=1.5) # The distribution of sd for supplementary data sources
sup_N <- c(15:25) # The sample size range for supplementary data sources
sup_homo <- 0 # The homogeneous (exchangeable) supplementry mean
sup_hetero_mean <- NULL

simulation_homo_cluster_hard <- function(mtd='bsl',perfect=F,k=10,q=10){ 
  
  ###################################
  # mtd: method, including naive kmeans ('kmeans'), iMEM ('bsl'), and clustering methods ('topk', 'even', 'split')
  # perfect: use changepoint detector or not
  # k: number of clusters for kmeans
  # q: number of selected sources for iMEM
  ###################################  
  
  # results are generated from 1000 random simulations
  for (i in 1:1000){
    
    num_homo <- 60
    
    # sampling primary source
    primary <- true_data[sample(1:10000, 20)] # The primary source is sampled from the true data with sample size = 20
    
    # sampling supplementary sources
    supple_gen_mean <- c(rep(sup_homo,num_homo),rep(0.5,10),rep(1,10),rep(1.5,10),rep(2,10))
    supple_gen_sd <- sup_sd[sample(1:10000,100)]
    supple_gen_N <- sup_N[sample(1:length(sup_N),100,replace = T)]
    
    # results for kmeans, bsl (iMEM), or dMEM
    if (mtd=='kmeans'){
      mem_result <- kmeans_baseline(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,perfect=perfect,k=k)
    } else if(mtd=='bsl'){
      mem_result <- clustering_methods(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,mtd=mtd,perfect=perfect,max_num_cluster=q)
    } else {
      mem_result <- clustering_methods(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,mtd=mtd,perfect=perfect)
    }
    
    # generate result table
    if (i==1){
      df <- t(c(mem_result['postmean'],mem_result['postvar'],mem_result['ess']))
      df <- as.data.frame(df)
      colnames(df) <- c('postmean','postvar','postess')
      
      write.csv(df,paste0('optimal_clustering_simulation_hard_',mtd,'method_',perfect,'perfect_',min(k,q),'kq.csv'),row.names=F)
      
    } else {
      df <- t(c(mem_result['postmean'],mem_result['postvar'],mem_result['ess']))
      write.table(df,paste0('optimal_clustering_simulation_hard_',mtd,'method_',perfect,'perfect_',min(k,q),'kq.csv'),append=T,sep=',',row.names=F,col.names=F)
    }
  }
  
}

simulation_random_cluster_hard  <- function(m,perfect=F, max_num_cluster=10){ 
  
  ###################################
  # m: times of randomly clustering within one iteration, results are the average of the m times
  # perfect: use changepoint detector or not
  # max_num_cluster: number of clusters
  ###################################
  
  # results are generated from 1000 random simulations
  for (i in 1:1000){
    
    num_homo <- 60
    
    # sampling primary source
    primary <- true_data[sample(1:10000, 20)] # The primary source is sampled from the true data with sample size = 20
    
    # sampling supplementary sources
    supple_gen_mean <- c(rep(sup_homo,num_homo),rep(0.5,10),rep(1,10),rep(1.5,10),rep(2,10))
    supple_gen_sd <- sup_sd[sample(1:10000,100)]
    supple_gen_N <- sup_N[sample(1:length(sup_N),100,replace = T)]
    
    # dMEM results for random clustering
    mem_result <- random_clustering(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,m=m,perfect=perfect, max_num_cluster=max_num_cluster)
    
    # generate result table
    if (i==1){
      df <- t(c(mem_result['postmean'],mem_result['postvar'],mem_result['ess']))
      df <- as.data.frame(df)
      colnames(df) <- c('postmean','postvar','postess')
      
      write.csv(df,paste0('optimal_clustering_simulation_hard_',m,'iter_',perfect,'perfect_',max_num_cluster,'max_clst.csv'),row.names=F)
      
    } else {
      df <- t(c(mem_result['postmean'],mem_result['postvar'],mem_result['ess']))
      write.table(df,paste0('optimal_clustering_simulation_hard_',m,'iter_',perfect,'perfect_',max_num_cluster,'max_clst.csv'),append=T,sep=',',row.names=F,col.names=F)
    }
  }
  
}


simulation_homo_cluster_hard(mtd='kmeans')
simulation_homo_cluster_hard(mtd='kmeans',k=5)
simulation_homo_cluster_hard(mtd='bsl')
simulation_homo_cluster_hard(mtd='bsl',q=5)
simulation_homo_cluster_hard(mtd='topk')
simulation_homo_cluster_hard(mtd='even')
simulation_homo_cluster_hard(mtd='split')
simulation_random_cluster_hard(1)
simulation_random_cluster_hard(10)
simulation_random_cluster_hard(1,perfect=F, max_num_cluster=1)
simulation_random_cluster_hard(1,perfect=F, max_num_cluster=3)


############################Simulation Scenaio 3 ############################
############### Varying proportions of nonexchangeable sources ##############
#############################################################################

# primary source
true_data <- rnorm(10000,0,sd=1) # The true distribution N(0,sd=1) for primary data source

# supplementary sources
sup_homo <- 0
sup_hetero_mean <- 1
sup_sd <- runif(50000,min=0.5,max=1.5) # The distribution of sd for supplementary data sources
sup_N <- c(15:25) # The sample size for supplementary data sources


simulation_proportion <- function(prop,num_supple,mtd='dMEM'){
  
  ###################################
  # prop: proportion of nonexchangeable sources
  # num_supple: total number of supplementary sources
  # mtd: iMEM (q=10) or dMEM (random averaging with m = 10)
  ###################################
  
  # results are generated from 1000 random simulations
  for (i in 1:1000){
    
    # sampling primary source
    primary <- true_data[sample(1:length(true_data), 20)] # The primary source is sampled from the true data with sample size = 20
    
    # sampling supplementary source
    supple_gen_mean <- c(rep(sup_homo,ceiling(prop*num_supple)),rep(sup_hetero_mean,num_supple-ceiling(prop*num_supple)))
    supple_gen_sd <- sup_sd[sample(1:length(sup_sd),num_supple,replace = T)]
    supple_gen_N <- sup_N[sample(1:length(sup_N),num_supple,replace = T)]
    
    # generate result & table
    if (mtd=='dMEM'){
      
      # dMEM with random averaging
      dmem_result <- random_clustering(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,m=10,perfect=F,max_num_cluster=10)
      
      # generate result table                           
      if (i==1){
        df <- t(c(dmem_result['postmean'],dmem_result['postvar'],dmem_result['ess']))
        df <- as.data.frame(df)
        colnames(df) <- c('postmean','postvar','postess')
        
        write.csv(df,paste0('simulation_homo_proprotion_',prop,'.csv'),row.names=F)
        
      } else {
        df <- t(c(dmem_result['postmean'],dmem_result['postvar'],dmem_result['ess']))
        write.table(df,paste0('simulation_homo_proprotion_',prop,'.csv'),append=T,sep=',',row.names=F,col.names=F)
      }
      
    } else if (mtd=='iMEM'){
      
      # iMEM with q = 10
      imem_result <- clustering_methods(primary,supple_gen_mean,supple_gen_sd,supple_gen_N,mtd='bsl',perfect=F)
      
      # generate result table
      if (i==1){
        df <- t(c(imem_result['postmean'],imem_result['postvar'],imem_result['ess']))
        df <- as.data.frame(df)
        colnames(df) <- c('postmean','postvar','postess')
        
        write.csv(df,paste0('simulation_homo_proprotion_',prop,'_imem.csv'),row.names=F)
        
      } else {
        df <- t(c(imem_result['postmean'],imem_result['postvar'],imem_result['ess']))
        write.table(df,paste0('simulation_homo_proprotion_',prop,'_imem.csv'),append=T,sep=',',row.names=F,col.names=F)
      }
    }
    
  }
  
}

lapply(seq(0.1,0.5,0.1),function(x){simulation_proportion(prop=x,num_supple=500)})
lapply(seq(0.1,0.5,0.1),function(x){simulation_proportion(prop=x,num_supple=500,mtd='iMEM')})

