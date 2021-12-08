source('iMEM_functions.R')

library(ggplot2)
library(scales)
library(gtable)
library(gridExtra)
library(ggtext)

simulation_optimal <- function(n,p,q,bsl=F,homo=60){ 
  #############################################################
  # n=number of clusters
  # p=length n vector of # homo grps in each cluster
  # q=length n vector of # hetero grps in each cluster
  # bsl=T/F, baseline
  # homo=number of truly homogeneous grps in simulated data
  ##############################################################
  
  # invalid p & q choice
  if(length(p)!=n | length(q)!=n){
    break
  }
  if(min(p+q)==0){
    break
  }
  if(sum(p)>homo|sum(q)>(100-homo)){
    break
  }
  
  for (i in 1:1000){
    
    # sample the primary and supplementary sources
    primary <- true_data[sample(1:10000, 20)] # The primary source is sampled from the true data with sample size = 10
    supple_gen_mean <- c(rep(sup_homo,homo),sup_hetero_mean[sample(1:10000,(100-homo))])
    supple_gen_sd <- sup_sd[sample(1:10000,100)]
    supple_N <- sup_N[sample(1:length(sup_N),100,replace = T)]
    
    supple <- list()
    for (k in 1:100){
      supple[[k]] <- c(rnorm(supple_N[k],supple_gen_mean[k],supple_gen_sd[k]))
    }
    
    if(bsl==T){ 
      # iMEM top 10
      supple_mean <- unlist(lapply(supple,mean))
      supple_sd <- unlist(lapply(supple,sd))
      mem_result <- imem_marg(c(mean(primary),sd(primary),length(primary)),
                              supple_mean[1:homo],supple_sd[1:homo],supple_N[1:homo],
                              prior='pi_e',final_grpsize = 10)

    } else {
      # MEM results with clustering
      supple_mean_cluster <- c()
      supple_sd_cluster <- c()
      supple_N_cluster <- c()
      
      for (j in 1:n){
        if(p[j]==0 & q[j]>0){
          # the cluster has only heterogeneous sources
          supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[(sum(q[0:(j-1)])+homo+1):(sum(q[1:j])+homo)])))
          supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[(sum(q[0:(j-1)])+homo+1):(sum(q[1:j])+homo)])))
          supple_N_cluster <- c(supple_N_cluster, sum(supple_N[(sum(q[0:(j-1)])+homo+1):(sum(q[1:j])+homo)]))
        } else if (q[j]==0 & p[j]>0){
          # the cluster has only homogeneous sources
          supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[(sum(p[0:(j-1)])+1):sum(p[1:j])])))
          supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[(sum(p[0:(j-1)])+1):sum(p[1:j])])))
          supple_N_cluster <- c(supple_N_cluster, sum(supple_N[(sum(p[0:(j-1)])+1):sum(p[1:j])]))
        } else {
          # the cluster is a mix of both homogeneous and heterogeneous sources
          supple_mean_cluster <- c(supple_mean_cluster,mean(unlist(supple[c((sum(p[0:(j-1)])+1):sum(p[1:j]),
                                                                            (sum(q[0:(j-1)])+homo+1):(sum(q[1:j])+homo))])))
          supple_sd_cluster <- c(supple_sd_cluster,sd(unlist(supple[c((sum(p[0:(j-1)])+1):sum(p[1:j]),
                                                                      (sum(q[0:(j-1)])+homo+1):(sum(q[1:j])+homo))])))
          supple_N_cluster <- c(supple_N_cluster, sum(supple_N[c((sum(p[0:(j-1)])+1):sum(p[1:j]),
                                                                 (sum(q[0:(j-1)])+homo+1):(sum(q[1:j])+homo))]))
        }
        
      }
      
      # MEM result
      mem_result <- mem_calc(c(mean(primary),sd(primary),20),
                             supple_mean_cluster,
                             supple_sd_cluster,
                             supple_N_cluster,
                             prior='pi_e')
      
    }
    
    # generate result table
    if (i==1){
      df <- t(c(mem_result$postmean,mem_result$postvar,mem_result$ess))
      df <- as.data.frame(df)
      colnames(df) <- c('postmean','postvar','postess')
      
      write.csv(df,paste0('optimal_simulation_',homo,'numhomo_',mean_hetero_mean,'heteromean_',n,'cluster_',paste0(p,collapse=""),'homoclst_',paste0(q,collapse=""),'heteroclst.csv'),row.names=F)
      
    } else {
      df <- t(c(mem_result$postmean,mem_result$postvar,mem_result$ess))
      write.table(df,paste0('optimal_simulation_',homo,'numhomo_',mean_hetero_mean,'heteromean_',n,'cluster_',paste0(p,collapse=""),'homoclst_',paste0(q,collapse=""),'heteroclst.csv'),append=T,sep=',',row.names=F,col.names=F)
    }
  }
  
}


mean_hetero_mean <- 1
true_data <- rnorm(10000,0,sd=1) # The true distribution N(0,sd=1) for primary data source
sup_sd <- rep(1,10000) # The SD for supplementary data sources
sup_N <- 20 # The sample size for supplementary data sources
sup_homo <- 0  # The homogeneous mean
sup_hetero_mean <- rep(mean_hetero_mean,10000)  # The heterogeneous mean


# running the simulations
simulation_optimal(n=10,p=c(rep(1,10)),q=c(rep(0,10)),bsl=T)
simulation_optimal(n=10,p=c(rep(10,6),rep(0,4)),q=c(rep(0,6),rep(10,4)))
simulation_optimal(n=10,p=c(rep(6,10)),q=c(rep(4,10)))
simulation_optimal(n=10,p=c(rep(6,10)),q=c(rep(0,10)))


library(plyr)
filename <- list.files(path = "../data/Motivating example results/",pattern = '.csv')
detach(package:plyr)

postvar <- c()
postvar_upper<-c()
postvar_lower<-c()
postvar_75<-c()
postvar_25<-c()
MSE <-c()
bias <- c()
MSE_75 <- c()
MSE_25 <- c()
MSE_upper <- c()
MSE_lower <- c()
bias_75 <- c()
bias_25 <- c()
bias_upper <- c()
bias_lower <- c()
ESS <- c()
ESS_75 <- c()
ESS_25 <- c()
ESS_upper <- c()
ESS_lower <- c()


for(i in 1:length(filename)) {
  data <- read.csv(paste0("../data/Motivating example results/",filename)[i])
  
  postvar <- c(postvar,median(data$postvar))
  postvar_upper<-c(postvar_upper,quantile(data$postvar,0.975))
  postvar_lower<-c(postvar_lower,quantile(data$postvar,0.025))
  postvar_75<-c(postvar_75,quantile(data$postvar,0.75))
  postvar_25<-c(postvar_25,quantile(data$postvar,0.25))
  
  MSE_upper<- c(MSE_upper,quantile((data$postmean)^2+data$postvar,0.975))
  MSE_lower<-c(MSE_lower,quantile((data$postmean)^2+data$postvar,0.025))
  MSE_75<- c(MSE_75,quantile((data$postmean)^2+data$postvar,0.75))
  MSE_25<-c(MSE_25,quantile((data$postmean)^2+data$postvar,0.25))
  MSE <- c(MSE,median((data$postmean)^2+data$postvar,na.rm=T))
  
  bias <- c(bias,median(data$postmean))
  bias_75<- c(bias_75,quantile((data$postmean),0.75))
  bias_25<-c(bias_25,quantile((data$postmean),0.25))
  bias_upper<- c(bias_upper,quantile((data$postmean),0.975))
  bias_lower<-c(bias_lower,quantile((data$postmean),0.025))
  
  ESS <- c(ESS,median(data$postess))
  ESS_75<- c(ESS_75,quantile((data$postess),0.75))
  ESS_25<-c(ESS_25,quantile((data$postess),0.25))
  ESS_upper<- c(ESS_upper,quantile((data$postess),0.975))
  ESS_lower<-c(ESS_lower,quantile((data$postess),0.025))
  
}

postvar_diff <- postvar-postvar[2]
case <- c('10x10','iMEM','10x6','10x(6+4)')
x<-c(3.5,1,4.75,2.25)

df <- as.data.frame(cbind(postvar_diff,postvar,postvar_upper,postvar_lower,postvar_75,postvar_25,MSE,MSE_upper,MSE_lower,MSE_75,MSE_25,bias,bias_upper,bias_lower,bias_75,bias_25,ESS,ESS_upper,ESS_lower,ESS_75,ESS_25,x))

pMSE <- ggplot(data=df,aes(x=x))+
  geom_boxplot(aes(ymin=MSE_lower,ymax=MSE_upper,lower=MSE_25,upper=MSE_75,middle=MSE,group=x),
               stat = "identity")+
  scale_x_continuous(breaks=c(3.5,1,4.75,2.25),labels=c('10x10','iMEM','10x6','10x(6+4)')) +
  labs(y='MSE',x='Clustering methods')+theme(axis.text.x = element_text(angle=30,hjust=1))

pESS <- ggplot(data=df,aes(x=x))+
  geom_boxplot(aes(ymin=ESS_lower,ymax=ESS_upper,lower=ESS_25,upper=ESS_75,middle=ESS,group=x),
               stat = "identity")+
  scale_x_continuous(breaks=c(3.5,1,4.75,2.25),labels=c('10x10','iMEM','10x6','10x(6+4)')) +
  labs(y='ESSS',x='Clustering methods')+theme(axis.text.x = element_text(angle=30, hjust=1))

grid.arrange(pMSE, pESS, ncol =2)

