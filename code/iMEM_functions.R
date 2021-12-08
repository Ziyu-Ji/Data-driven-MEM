
# This file is the open source implementation of iMEM by Roland Brown et al
# The GitHub source is https://github.com/rbrown789/Iterated-Multisource-Exchangeability-Models



############## FUNCTIONS FOR IMPLEMENTING STANDARD GAUSSIAN MEMS ########################## 

modpostmean <- function(modvec,prim, means, Vs)
{
  ########################################
  # calculate posterior mean for a single exchangeability model specified by modvec
  #
  #  modvec: logical vector specifying a particular exchangeability model
  #  prim: infomration on the primary source (format: c(mean,sd,N) )
  #  means: means for all supplementary sources (must be in same order as modvec)
  #  Vs: s^2/n for all supplmentary sources (must be in same order as modvec)
  #	
  ########################################
  
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  # calculate posterior mean from Eq. 3.6 in paper
  pmean <- (m*prod(Vs^modvec) + sum( (v*means*modvec)/Vs )*prod(Vs^modvec) ) / 
    ( v*( sum(modvec/Vs)*prod(Vs^modvec) ) + prod(Vs^modvec) ) 
  
  return(pmean)
}

modpostvar <- function(modvec,prim, means, Vs)
{
  ########################################
  # calculate posterior variance for a single exchangeability model specified by modvec
  #
  #  modvec: logical vector specifying a particular exchangeability model
  #  prim: infomration on the primary source (format: c(mean,sd,N) )
  #  means: means for all supplementary sources (must be in same order as modvec)
  #  Vs: s^2/n for all supplmentary sources (must be in same order as modvec)
  #	
  ########################################
  
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  # calculate posterior variance from Eq. 3.6 in paper
  pvar <- ( 1/v + sum( modvec/Vs ) )^-1	
  
  return( pvar )	
}



modmlik <- function(modvec,prim,means,Vs)
{
  ################################################################
  # calculate conditional marginal likelhiood for single exchangeability model specified by modvec 
  #
  #  modvec: logical vector specifying a particular exchangeability model
  #  prim: infomration on the primary source (format: c(mean,sd,N) )
  #  means: means for all supplementary sources (must be in same order as modvec)
  #  Vs: s^2/n for all supplmentary sources (must be in same order as modvec)
  #
  ################################################################
  
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  # number of total sources
  nsrc <- length(means)
  
  #########
  
  # calculate the marginal likelihood given the exchangebility model specified by modvec (per Eq. 3.5 in paper)
  
  # calculate internal sum A needed for marginal likelihood computation
  As <- rep(NA,length(modvec))
  for( i in 1:length(modvec)) { As[i] <- sum( modvec[-i]/Vs[-i] ) }
  
  # calculate internal sum B needed for marginal likelihood computation
  Bs <- rep(NA,length(modvec))
  for(i in 1:length(modvec)) 
  {
    gtind <- (1:length(modvec)) > i
    
    modvec2 <- modvec[gtind]
    means2 <- means[gtind]
    Vs2 <- Vs[gtind]
    
    # calculate values for the sum within quantity B
    Cs <- rep(NA,length( (i+1):length(modvec) ))
    k <- 1
    for(j in (i+1):length(modvec)) 
    {
      Cs[k] <- sum( (modvec/Vs)[-c(i,j)] )
      k <- k+1
    }		
    
    Bs[i] <- sum( (modvec[i]*modvec2*(means[i]-means2)^2)/(Vs[i] + Vs2 + Vs[i]*Vs2*( 1/v + Cs) ))		
  }	
  
  mlik <- ( sqrt(2*pi)^(nsrc + 1 - sum(modvec) ) / sqrt( (1/v + sum(modvec/Vs))*prod( (1/Vs)^!modvec )) )*
    exp( -0.5*sum( (modvec*(m-means)^2)/(v+Vs+v*Vs*As ) + Bs ) )		
  
  return(mlik)	
}



modweights <- function(mmat,mliks,prim,means,Vs,prior)
{
  
  ################################################################
  # calculate weights for all exchangeability models specified in mmat 
  #
  #  mmat: logical matrix of all possible exchangeability models (rows: models, columns: sources)
  #  mliks: vector of conditional marginal likelhioods for each possible model (must be in same order as cols of mmat)
  #  prim: infomration on the primary source (format: c(mean,sd,N) )
  #  means: means for all supplementary sources (must be in same order as rows of mmat)
  #  Vs: s^2/n for all supplmentary sources (must be in same order as rows of mmat)
  #
  ################################################################
  
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  # number of total sources
  nsrc <- length(means)
  
  #########
  
  # equal prior weights for pi_e
  if( prior=='pi_e') { modpriors <- 1 
  
  # for pi_n use Eq. 3.7 from paper
  }else if (prior=='pi_n') {
    
    # quantities needed for the source prior calcualtion
    As <- apply(mmat,1, function(modvec) { sum( modvec/Vs) } )
    Bs <- apply(mmat,1, function(modvec) { prod( (1/Vs)^(!modvec) ) } )
    Cs <- apply(mmat,1, function(modvec) { sum( modvec ) } )
    
    ### need to check if it's 1/s or 1/s^2, either typo in paper, 
    ### or mistake in Alex's code, currently set to match Alex's code
    Ds <- sqrt( (1/s + As)*Bs)/ ( sqrt(2*pi)^( nsrc + 1- Cs )) 
    
    # calculate source inclusion priors
    priors_1 <- apply( mmat,2,function(srcvec) { sum( srcvec*Ds )  } )
    priors_0 <- apply( mmat,2,function(srcvec) { sum( (!srcvec)*Ds) } )
    
    # calculate model inclusion priors from source inclusion priors
    modpriors <- apply(mmat,1,function(modvec){ prod(modvec*priors_1 + (!modvec)*priors_0) })
    
  } else { stop("Current support only for priors 'pi_e' and 'pi_n") }
  
  # calculate raw posterior weights
  rawwts <- modpriors*mliks
  
  # calculate scaled posterior weights
  wts <- rawwts/sum(rawwts)
  
  return(wts) 	
}


derivexp <- function(modvec,prim, means, Vs)
{
  ########################################
  # calculate derivative of exponential term in the marginal likelhiood
  #  for single exchangeability model specified by modvec
  #
  #  modvec: logical vector specifying a particular exchangeability model
  #  prim: infomration on the primary source (format: c(mean,sd,N) )
  #  means: means for all supplementary sources (must be in same order as modvec)
  #  Vs: s^2/n for all supplmentary sources (must be in same order as modvec)
  #	
  ########################################
  
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  # calculate internal sum A 
  As <- rep(NA,length(modvec))
  for( i in 1:length(modvec)) { As[i] <- sum( modvec[-i]/Vs[-i] ) }
  
  # calcualte
  derexp <- sum ( (-modvec*(m-means))/(v+Vs+v*Vs*As ) )
  return( derexp )	
}



derivwts <- function(mods)
{
  ########################################
  # calculate derivative of weights for single exchangeability model specified by modvec
  #
  #  mods: the 
  #	
  ########################################
  
  num <- rep(NA,nrow(mods))
  for(i in 1:length(num) ) { num[i] <- sum( mods$postwts*(mods$expder- mods$expder[i]) )}
  return(num/mods$postwts)
  
}




derivpmean <- function(modvec,prim,Vs)
{
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  return( prod(Vs^modvec)/( v*( sum(modvec/Vs)*prod(Vs^modvec) ) + prod(Vs^modvec) ) )
}





############################################################



mem_calc <- function( prim, means,sds, Ns, prior='pi_e' )
{
  ##########################################################################
  # This function carries out MEM calcualtions for the case of gaussian means for an arbitrary number
  # of supplmentary sources.  All calculations are closed form, but model space increases exponentially
  # with linearly increasing number of sources, so unsure how efficiently it will scale 
  # 
  #   prim: infomration on the primary source (format: c(mean,sd,N) )
  #   means: vector of means for supplementary sources
  #   sds: vector of standard deviations for supplementary sources
  #   Ns: vector of sample sizes for supplementary sources
  #   prior: type of prior
  #  
  # NOTE: means, sds, and Ns must all be the same length
  ############################################################################
  
  # if means, sds, and Ns are of different length, error out
  if( length(unique( c(length(means),length(sds),length(Ns)) )) > 1 ){ stop("'means','sds','Ns' must be of equal length")}
  
  # First, let's generate a matrix with all possible exchangebility models
  nsrc <- length(means)	
  mods <- expand.grid ( rep(list(c(FALSE,TRUE)), nsrc)); names(mods) <- paste0("src",1:nsrc) 
  mmat <- as.matrix ( mods )
  
  # calculate v = s^2/n for suppelementary sources
  Vs <- sds^2/Ns
  
  # Calculate model-specific posterior means
  mods$pmeans <- apply(mmat,1,modpostmean, prim, means,Vs)
  
  # Calculate model-specific posterior variances
  mods$pvars <- apply(mmat,1,modpostvar, prim, means,Vs)
  
  # calculate conditional marginal likelihood for each model
  mods$mlik <- apply(mmat,1,modmlik,prim,means,Vs)
  
  # calculate posterior weights for each model
  mods$postwts <- modweights(mmat,mods$mlik,prim,means,Vs,prior)
  
  ## calculate overall posterior mean ##
  pmean <- sum(mods$pmeans*mods$postwts)
  
  ### calcualte overall posterior variance ###
  mods$expder <- apply(mmat,1,derivexp,prim,means,Vs) # derivative of exponential part of marginal likelhiood
  mods$wtsder <- derivwts(mods) # derivative of weight calculation
  mods$pmeander <- apply(mmat,1,derivpmean,prim,Vs) # derivative of posterior mean
  mods$gder <- ( (mods$pmeander/mods$postwts) - (mods$pmeans*mods$wtsder) )/ (1/mods$postwts)^2
  mods$gder[is.nan(mods$gder)] <- 0
  gmat <- mods$gder%*%t(mods$gder)
  pvar <- rep(1,nrow(mods))%*%gmat%*%rep(1,nrow(mods))* ( (prim[2]^2)/prim[3])
  
  ### calculate effective historical sample sizes for individual models, and overall
  mods$ess <- (prim[2]^2/mods$pvars)-prim[3]
  ess_fin <- sum(mods$ess*mods$postwts)
  
  
  #############################
  
  # generate output list and return
  srcinfo <- cbind( prim, rbind(means,sds,Ns) )
  colnames(srcinfo) <- c("primary",paste0("src",1:(ncol(srcinfo)-1)))
  rownames(srcinfo) <- c("mean","sd","N")
  
  out <- list(postmean=pmean,postvar=pvar,ess=ess_fin,
              memlist=mods[,c(grep("src",names(mods),value=T),"pmeans","pvars","ess","postwts")],
              srcinfo=srcinfo)
  return(out)
}

#######################################################################################################
#######################################################################################################



####### MARGINAL IMEM FUNCTIONS #######

imem_marg <- function(prim, means,sds, Ns, prior='pi_e',final_grpsize,qthresh)
{
  ####################################################################
  # function to fit marginal selection iterated MEM 
  #
  # prim: primary source information c(mean,sd,N)
  # means: supplementary  means
  # sds: supplementary sds
  # Ns: supplementary sample sizes
  # prior: prior desired, 'pi_e' or 'pi_n'
  # final_grpsize: q, the number of sources in the final MEM model, must be <= # supp. sources
  #
  #
  ###################################################################
  
  
  # final_grpsize <- 3
  # prim <- c(-2,4,20)
  # means <- xbars3
  # sds <- rep(sigma,length(means))
  # Ns <- rep(n,length(means))
  
  # first fit all marginal MEM models to be used in mem_calc()
  margmems <- lapply(1:length(means),function(x) { mem_calc(prim=prim,means=means[x],sds=sds[x],Ns=Ns[x],prior='pi_e' )} )
  
  # now extract the scores and order them
  margscores <- sapply(1:length(margmems), function(x){ margmems[[x]]$memlist$postwts[2] })
  ord <- rev(order(margscores))
  
  # if qthresh is specified
  if(final_grpsize=="thresh") {
    final_grpsize <- sum(margscores >= qthresh)
    if(final_grpsize > 12) { final_grpsize <- 12}
  }
  
  
  # print(final_grpsize)
  
  # grab the sources with q highest scores
  finmeans <- means[ord][1:final_grpsize]
  finsds <- sds[ord][1:final_grpsize]
  finNs <- Ns[ord][1:final_grpsize]
  finScores <- margscores[ord][1:final_grpsize]
  
  # fit the final mem model
  finmem <- mem_calc(prim=prim,means=finmeans,sds=finsds,Ns=finNs,prior=prior)
  finmem$srcinfo <- rbind(finmem$srcinfo,c(NA,finScores) )
  rownames(finmem$srcinfo)[4] <- "score"
  
  # if final_grpsize is zero
  if(final_grpsize==0){
    
    finmem$postmean <- finmem$memlist$pmeans[1]
    finmem$postvar <- finmem$memlist$pvars[1]
    finmem$ess <- 0
    finmem$memlist <- finmem$memlist[1,c("pmeans","pvars","ess","postwts")]
    finmem$memlist$postwts <- 1
    finmem$srcinfo <- as.matrix( finmem$srcinfo[,1],colnames="primary")
    colnames(finmem$srcinfo) <- "primary"
  }
  
  return(finmem)
}


imem_marg_v2 <- function(prim, means,sds, Ns, IDs,prior='pi_e',final_grpsize)
{
  ####################################################################
  # Same functionality as imem_marg(), but additionally retains the IDs of the sources
  #  selected into the final model.  Additional input variable is `IDs`, which are 
  #  identifiers for the supplementary sources
  #
  ###################################################################
  
  
  # final_grpsize <- 3
  # prim <- c(-2,4,20)
  # means <- xbars3
  # sds <- rep(sigma,length(means))
  # Ns <- rep(n,length(means))
  
  # first fit all marginal MEM models to be used in mem_calc()
  margmems <- lapply(1:length(means),function(x) { mem_calc(prim=prim,means=means[x],sds=sds[x],Ns=Ns[x],prior='pi_e' )} )
  
  # now extract the scores and order them
  margscores <- sapply(1:length(margmems), function(x){ margmems[[x]]$memlist$postwts[2] })
  ord <- rev(order(margscores))
  
  # grab the highest "final_grpsize" highest scores
  finmeans <- means[ord][1:final_grpsize]
  finsds <- sds[ord][1:final_grpsize]
  finNs <- Ns[ord][1:final_grpsize]
  finIDs <- IDs[ord][1:final_grpsize]
  finScores <- margscores[ord][1:final_grpsize]
  
  # fit the final mem model
  finmem <- mem_calc(prim=prim,means=finmeans,sds=finsds,Ns=finNs,prior=prior)
  finmem$ids.final <- finIDs
  finmem$srcinfo <- rbind(finmem$srcinfo,c(NA,finScores) )
  rownames(finmem$srcinfo)[4] <- "score"
  colnames(finmem$srcinfo) <- c("primary",finIDs)
  colnames(finmem$memlist)[1:(ncol(finmem$memlist)-4)] <- finIDs
  
  return(finmem)
}


