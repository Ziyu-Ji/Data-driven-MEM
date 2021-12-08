source('iMEM_functions.R')
source('dMEM_functions.R')

# Function to grab all trips/activities of a certain type for a certain i
typeid <- function(dat,id,type,subtype,varnms=names(dat))
{
  return(dat[dat$surveycode %in% id & dat$type %in% type & dat$primary_mode %in% subtype, varnms])
}


# iMEM function given subtype, primary ID, and dataset for supplementary ids with > 5 observations
imem_dayn <- function(outname,type,subtype,primid,nfinsrc,qthresh=NULL,dat)
{
  # create primary source information vector
  prim_dat <- typeid(dat,primid,type,subtype,varnms=c("surveycode",outname))
  prim_out <- prim_dat[[outname]]
  prim_info <- c(mean(prim_out,na.rm=T),sd(prim_out,na.rm=T),sum(!is.na(prim_out))  )
  
  if(prim_info[2] ==0) { 
    out <- list(postmean=prim_info[1],postvar=0,note="standard deviation 0 for primary source")
  } else {
    # find all the valid supplementary ids with >= 5 observations
    tmp <- typeid(dat,unique(dat$surveycode)[!unique(dat$surveycode) %in% primid],type,subtype,varnms=c("surveycode",outname))
    tmp <- tmp[!is.na(tmp[[outname]]),]
    tbl <- rev(sort(table(tmp$surveycode)))
    suppids <- names(tbl)[tbl >= 5]
    
    # create supplmentary information vectors
    suppmeans <- suppsds <- suppNs <- suppids2 <- c()
    suppdat_all <- NULL
    
    for(i in suppids){
      supp_dat <- typeid(dat,i,type,subtype,varnms=c("surveycode",outname) )
      supp_out <- supp_dat[[outname]]
      
      tmp.sd <- sd(supp_out,na.rm=T)
      
      if(tmp.sd > 0) {
        suppmeans <- c(suppmeans,mean(supp_out,na.rm=T) )
        suppsds <- c(suppsds,tmp.sd)
        suppNs <- c(suppNs,sum(!is.na(supp_out)))
        suppids2 <- c(suppids2,i)
        suppdat_all <- rbind(suppdat_all,supp_dat)
      }
    }
    
    # apply MEM function to the primary and supplmentary information
    out <- imem_marg(prim=prim_info,means=suppmeans,sds=suppsds,Ns=suppNs,prior="pi_e",final_grpsize = nfinsrc,qthresh=qthresh)
  }
  
  return(out)
}

# dMEM function given subtype, primary ID, and dataset for supplementary ids with > 5 observations
DMEM_dayn <- function(outname,type="TRIP",subtype,primid,mtd,m=10,max_num_cluster=10,dat,lst=F){
  # create primary source information vector
  prim_dat <- typeid(dat,primid,type,subtype,varnms=c("surveycode",outname))
  prim_out <- prim_dat[[outname]]
  prim_info <- c(mean(prim_out,na.rm=T),sd(prim_out,na.rm=T),sum(!is.na(prim_out))  )
  
  # find all the valid supplementary ids with >= 5 observations
  tmp <- typeid(dat,unique(dat$surveycode)[!unique(dat$surveycode) %in% primid],type,subtype,varnms=c("surveycode",outname))
  tmp <- tmp[!is.na(tmp[[outname]]),]
  tbl <- rev(sort(table(tmp$surveycode)))
  suppids <- names(tbl)[tbl >= 5]
  
  # create supplmentary information vectors
  suppmeans <- suppsds <- suppNs <- suppids2 <- c()
  suppdat_all <- NULL
  
  for(i in suppids){
    supp_dat <- list(na.omit(typeid(dat,i,type,subtype,varnms=c("surveycode",outname))[,2]))
    supp_out <- typeid(dat,i,type,subtype,varnms=c("surveycode",outname))[,2]
    
    tmp.sd <- sd(supp_out,na.rm=T)
    
    if(tmp.sd > 0) {
      suppmeans <- c(suppmeans,mean(supp_out,na.rm=T) )
      suppsds <- c(suppsds,tmp.sd)
      suppNs <- c(suppNs,sum(!is.na(supp_out)))
      suppids2 <- c(suppids2,i)
      suppdat_all <- append(suppdat_all,supp_dat)
    }
  }
  out <- DMEM_app(primary=prim_info,supple_mean=suppmeans,supple_sd=suppsds,supple_N=suppNs,supple=suppdat_all,
                  mtd=mtd,m=m,max_num_cluster=max_num_cluster)
  return(out)
}


# calculate posterior credible intervals from imem output
credint.imem <- function(imemobj,cred.lev=0.95)
{
  # sample the multinomial distribution for the normal mixture (out of 10000000 draws)
  nsamp.mod <- c(rmultinom(1,1000000,imemobj$memlist$postwts))
  
  # sample the normal mixture
  post.dist <- unlist(apply(cbind(nsamp.mod,imemobj$memlist$pmeans,imemobj$memlist$pvars), 1, function(x){rnorm(x[1],x[2],x[3])}))
  
  # calculate the HPD interval
  hpdint <- boa.hpd(post.dist,alpha=1-cred.lev)
  pm <- imemobj$postmean
  pse <- sqrt(imemobj$postvar)
  upr <- hpdint[1]
  lwr <- hpdint[2]
  esss <- imemobj$ess
  
  out <- c(pm,pse,lwr,upr,esss)
  names(out) <- c("imem.postmean","imem.postsd","imem.lwr","imem.upr","imem.esss")
  
  return(out)
}

# calculate posterior credible intervals from dmem output
credint.dmem <- function(dmemobj,cred.lev=0.95)
{
  dmemobj <- unlist(dmemobj)
  pm <- dmemobj[1]
  pse <- sqrt(dmemobj[2])
  esss <- dmemobj[3]
  upr <- dmemobj[4]
  lwr <- dmemobj[5]
  
  out <- c(pm,pse,lwr,upr,esss)
  names(out) <- c("dmem.postmean","dmem.postsd","dmem.lwr","dmem.upr","dmem.esss")
  
  return(out)
}

# calculate posterior credible intervals for simple mean result from imem output
credint.simplemean <- function(imemobj,cred.lev=0.95)
{
  pm <- imemobj$memlist$pmeans[1]
  pse <- sqrt(imemobj$memlist$pvars[1])
  upr <- pm + pse*qnorm( 1-(1-cred.lev)/2 )
  lwr <- pm - pse*qnorm( 1-(1-cred.lev)/2 )
  
  out <- c(pm,pse,lwr,upr)
  names(out) <- c("smean","ssd","slwr","supr")
  
  return(out)
}


# function to generate iMEM/dMEM/simple mean results for a set of primary ID
# select people with at least 2 different observations (to calculate sd) 
# primary sources are those have >= 10 observations in subtypes "WALK","CAR","BIKE","IN_VEHICLE" or "BUS"; 
# all supplementary sources are required to have >= 5 observations (built-in)
gendat.emotbymode <- function(primid,mtd,m=1,data)
{
  for(j in primid){
    # for each individual, drop the subtypes with <= 10 observations
    num <- data %>%
      filter(surveycode == j) %>%
      group_by(primary_mode) %>%
      summarise(num=length(type)) %>%
      filter(num >= 10)
    subtype <- intersect(num$primary_mode, c("WALK","CAR","BIKE","IN_VEHICLE","BUS"))
    
    if (length(subtype) > 0){
      for(k in subtype) {
        # for each subtype, supplementary source
        uniquedata <- data %>%
          filter(surveycode == j & primary_mode == k) %>%
          select(colnames(data)[-c(1:3)]) %>%
          gather(value="emo") %>%
          group_by(key) %>%
          summarise(sd=sd(emo,na.rm = T)) %>%
          filter(is.na(sd)==F)
        outcms <- uniquedata$key[uniquedata$sd!=0]
        
        if (length(outcms) > 0){
          for(i in outcms) {
            # fit the imem models
            imem.out <- imem_dayn(i,type="TRIP",subtype=k,primid=j,nfinsrc=10,dat=data)
            # generate plotting data with posterior means/intervals, simple means/intervals, esss
            plotdat <- data.frame(t(credint.imem(imem.out, cred.lev = 0.95)))
            
            # fit the dmem models
            dmem.out <- DMEM_dayn(i,type="TRIP",subtype=k,primid=j,mtd=mtd,m=m,max_num_cluster=10,dat=data)
            dmem.plotdat <- data.frame(t(credint.dmem(dmem.out, cred.lev=0.95)))
            plotdat <- cbind(plotdat,dmem.plotdat)
            
            # calculate the results for simple mean
            smean <-  data.frame(t(credint.simplemean(imem.out, cred.lev=0.95)))
            plotdat <- cbind(plotdat,smean)
            
            # format the output table
            plotdat$mode <- paste0(k,' (N=',round(imem.out$srcinfo[3,1],0),')',sep='')
            plotdat$outcm <- i
            plotdat$primid <- j
            plotdat$mode <- as.factor(plotdat$mode)
            plotdat$outcm <- as.factor(plotdat$outcm)
            plotdat$dmem.pctdiffmean <- abs(plotdat$dmem.postmean-plotdat$smean)/plotdat$smean
            plotdat$dmem.pctdecsd <- (plotdat$ssd-plotdat$dmem.postsd)/plotdat$ssd
            plotdat$imem.pctdiffmean <- abs(plotdat$imem.postmean-plotdat$smean)/plotdat$smean
            plotdat$imem.pctdecsd <- (plotdat$ssd-plotdat$imem.postsd)/plotdat$ssd
            
            if (j==primid[1] & k==subtype[1] & i==outcms[1]){
              write.csv(plotdat,paste0('application_',m,'m_',mtd,'method.csv'),row.names=F)
              
            } else {
              write.table(plotdat,paste0('application_',m,'m_',mtd,'method.csv'),append=T,sep=',',row.names=F,col.names=F)
            } 
          }
        }
      }
    }
  }
}


# function for QC check based on marginal weights
DMEM_qc <- function(primid,data)
{
  for(j in primid){
    num <- data %>%
      filter(surveycode == j) %>%
      group_by(primary_mode) %>%
      summarise(num=length(type)) %>%
      filter(num >= 10)
    subtype <- intersect(num$primary_mode, c("WALK","CAR","BIKE","IN_VEHICLE","BUS"))
    
    if (length(subtype) > 0){
      for(k in subtype) {
        uniquedata <- data %>%
          filter(surveycode == j & primary_mode == k) %>%
          select(colnames(data)[-c(1:3)]) %>%
          gather(value="emo") %>%
          group_by(key) %>%
          summarise(sd=sd(emo,na.rm = T)) %>%
          filter(is.na(sd)==F)
        outcms <- uniquedata$key[uniquedata$sd!=0]
        
        if (length(outcms) > 0){
          for(i in outcms) {
            qcdat <- NULL
            dmem <- DMEM_dayn(i,type="TRIP",subtype=k,primid=j,mtd=mtd,m=m,max_num_cluster=10,dat=data,result='qc')
            
            qcdat <- as.data.frame(rbind(dmem))
            
            qcdat$mode <- k
            qcdat$outcm <- i
            qcdat$primid <- j
            qcdat <- as.data.frame(qcdat)
            colnames(qcdat) <- c('selected supple','total supple','weights','mode','outcm','primid')
            qcdat <- t(apply(qcdat,2,as.character))
            if (j==primid[1] & i == outcms[1]){
              write.csv(qcdat,paste0('app_qc_all.csv'),row.names=F)
            } else {
              write.table(qcdat,paste0('app_qc_all.csv'),append=T,sep=',',row.names=F,col.names=F)
            }
          }
        }
      }
    }
  }
}
