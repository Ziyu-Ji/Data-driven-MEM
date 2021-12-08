library(ggplot2)
library(scales)
library(gtable)
library(gridExtra)
library(ggtext)

############################## Simulation plots#####################################

################################# scenario I #######################################
dir.path <- "../data/Simulation results/Scenario I/"
library(plyr)
filename <- list.files(path = "../data/Simulation results/Scenario I/",pattern = '.csv')
detach(package:plyr)

postvar <- c()
postvar_upper<-c()
postvar_lower<-c()
postvar_75<-c()
postvar_25<-c()
RMSE <-c()
bias <- c()
RMSE_75 <- c()
RMSE_25 <- c()
RMSE_upper <- c()
RMSE_lower <- c()
bias_75 <- c()
bias_25 <- c()
bias_upper <- c()
bias_lower <- c()
ESSS <- c()
ESSS_75 <- c()
ESSS_25 <- c()
ESSS_upper <- c()
ESSS_lower <- c()

for(i in 1:length(filename)) {
  data <- read.csv(paste0(dir.path,filename)[i])
  
  postvar <- c(postvar,median(data$postvar,na.rm=T))
  postvar_upper<-c(postvar_upper,quantile(data$postvar,0.975,na.rm=T))
  postvar_lower<-c(postvar_lower,quantile(data$postvar,0.025,na.rm=T))
  postvar_75<-c(postvar_75,quantile(data$postvar,0.75,na.rm=T))
  postvar_25<-c(postvar_25,quantile(data$postvar,0.25,na.rm=T))
  
  RMSE_upper<- c(RMSE_upper,quantile(sqrt((data$postmean)^2+data$postvar),0.975,na.rm=T))
  RMSE_lower<-c(RMSE_lower,quantile(sqrt((data$postmean)^2+data$postvar),0.025,na.rm=T))
  RMSE_75<- c(RMSE_75,quantile(sqrt((data$postmean)^2+data$postvar),0.75,na.rm=T))
  RMSE_25<-c(RMSE_25,quantile(sqrt((data$postmean)^2+data$postvar),0.25,na.rm=T))
  RMSE <- c(RMSE,median(sqrt((data$postmean)^2+data$postvar),na.rm=T))
  
  bias <- c(bias,median(data$postmean,na.rm=T))
  bias_75<- c(bias_75,quantile((data$postmean),0.75,na.rm=T))
  bias_25<-c(bias_25,quantile((data$postmean),0.25,na.rm=T))
  bias_upper<- c(bias_upper,quantile((data$postmean),0.975,na.rm=T))
  bias_lower<-c(bias_lower,quantile((data$postmean),0.025,na.rm=T))
  
  ESSS <- c(ESSS, mean(data$postess,na.rm=T))
  ESSS_75<- c(ESSS_75,quantile((data$postess),0.75,na.rm=T))
  ESSS_25<-c(ESSS_25,quantile((data$postess),0.25,na.rm=T))
  ESSS_upper<- c(ESSS_upper,quantile((data$postess),0.975,na.rm=T))
  ESSS_lower<-c(ESSS_lower,quantile((data$postess),0.025,na.rm=T))
  
}


postvar_diff <- postvar-postvar[1]
case <- c('random averaging','baseline: iMEM q=10','baseline: iMEM q=5', 'evenly distributed','baseline: kmeans k=10','baseline: kmeans k=5','top 5 + combine 5','combine top k', 'random','1 cluster','3 cluster')

order <- c('baseline: kmeans k=5','baseline: kmeans k=10','baseline: iMEM q=5','baseline: iMEM q=10','combine top k','evenly distributed','top 5 + combine 5', 'random', 'random averaging','3 cluster','1 cluster')

x <- c()
for (k in 1:length(case)){
  w <- which(order==case[k])
  plugin <- 1+(w-1)*1.25
  x <- c(x, plugin)
}

df1 <- as.data.frame(cbind(postvar_diff,postvar,postvar_upper,postvar_lower,postvar_75,postvar_25,RMSE,RMSE_upper,RMSE_lower,RMSE_75,RMSE_25,bias,bias_upper,bias_lower,bias_75,bias_25,ESSS,ESSS_upper,ESSS_lower,ESSS_75,ESSS_25,x))

s1var <- ggplot(data=df1,aes(x=x))+
  geom_boxplot(aes(ymin=postvar_lower,ymax=postvar_upper,lower=postvar_25,upper=postvar_75,middle=postvar,group=x),
               stat = 'identity', width = 0.8)+
  scale_x_continuous(breaks=seq(1,13.5,1.25),labels=c('kmeans k=5','kmeans k=10','iMEM q=5','iMEM q=10','ordered combine', 'evenly distributed', 'top half & combine half', 'random','random averaging','3 cluster','1 cluster')) +
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x))) +
  ylab('Posterior var')+xlab('(a) Scenario I')+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle=30,hjust=1)) 

s1bias <- ggplot(data=df1,aes(x=x))+
  geom_boxplot(aes(ymin=bias_lower,ymax=bias_upper,lower=bias_25,upper=bias_75,middle=bias,group=x),
               stat = 'identity', width = 0.8)+
  scale_x_continuous(breaks=seq(1,13.5,1.25),labels=c('kmeans k=5','kmeans k=10','iMEM q=5','iMEM q=10','ordered combine', 'evenly distributed', 'top half & combine half', 'random','random averaging','3 cluster','1 cluster')) + ylab('Bias') +xlab('(a) Scenario I')+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle=30,hjust=1))

s1RMSE <- ggplot(data=df1,aes(x=x))+
  geom_boxplot(aes(ymin=RMSE_lower,ymax=RMSE_upper,lower=RMSE_25,upper=RMSE_75,middle=RMSE,group=x),
               stat = 'identity', width = 0.8)+
  scale_x_continuous(breaks=seq(1,13.5,1.25),labels=c('kmeans k=5','kmeans k=10','iMEM q=5','iMEM q=10','ordered combine', 'evenly distributed', 'top half & combine half', 'random','random averaging','3 cluster','1 cluster')) +
  # scale_y_continuous(trans = 'log10',
  #                       breaks = trans_breaks('log10', function(x) 10^x),
  #                       labels = trans_format('log10', math_format(10^.x))) +
  ylab('RMSE')+xlab('(a) Scenario I')+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle=30,hjust=1))

s1esss <- ggplot(data=df1,aes(x=x))+
  geom_boxplot(aes(ymin=ESSS_lower,ymax=ESSS_upper,lower=ESSS_25,upper=ESSS_75,middle=ESSS,group=x),
               stat = 'identity', width = 0.8)+
  scale_x_continuous(breaks=seq(1,13.5,1.25),labels=c('kmeans k=5','kmeans k=10','iMEM q=5','iMEM q=10','ordered combine', 'evenly distributed', 'top half & combine half', 'random','random averaging','3 cluster','1 cluster')) + ylab('ESSS') +xlab('(a) Scenario I')+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle=30,hjust=1))

s1var
s1bias
s1RMSE
s1esss

################################# scenario II #######################################
dir.path <- "../data/Simulation results/Scenario II/"
library(plyr)
filename <- list.files(path = dir.path,pattern = '.csv')
detach(package:plyr)

postvar <- c()
postvar_upper<-c()
postvar_lower<-c()
postvar_75<-c()
postvar_25<-c()
RMSE <-c()
bias <- c()
RMSE_75 <- c()
RMSE_25 <- c()
RMSE_upper <- c()
RMSE_lower <- c()
bias_75 <- c()
bias_25 <- c()
bias_upper <- c()
bias_lower <- c()
ESSS <- c()
ESSS_75 <- c()
ESSS_25 <- c()
ESSS_upper <- c()
ESSS_lower <- c()


for(i in 1:length(filename)) {
  data <- read.csv(paste0(dir.path,filename)[i])
  
  postvar <- c(postvar,median(data$postvar,na.rm=T))
  postvar_upper<-c(postvar_upper,quantile(data$postvar,0.975,na.rm=T))
  postvar_lower<-c(postvar_lower,quantile(data$postvar,0.025,na.rm=T))
  postvar_75<-c(postvar_75,quantile(data$postvar,0.75,na.rm=T))
  postvar_25<-c(postvar_25,quantile(data$postvar,0.25,na.rm=T))
  
  RMSE_upper<- c(RMSE_upper,quantile(sqrt((data$postmean)^2+data$postvar),0.975,na.rm=T))
  RMSE_lower<-c(RMSE_lower,quantile(sqrt((data$postmean)^2+data$postvar),0.025,na.rm=T))
  RMSE_75<- c(RMSE_75,quantile(sqrt((data$postmean)^2+data$postvar),0.75,na.rm=T))
  RMSE_25<-c(RMSE_25,quantile(sqrt((data$postmean)^2+data$postvar),0.25,na.rm=T))
  RMSE <- c(RMSE,median(sqrt((data$postmean)^2+data$postvar),na.rm=T))
  
  bias <- c(bias,median(data$postmean,na.rm=T))
  bias_75<- c(bias_75,quantile((data$postmean),0.75,na.rm=T))
  bias_25<-c(bias_25,quantile((data$postmean),0.25,na.rm=T))
  bias_upper<- c(bias_upper,quantile((data$postmean),0.975,na.rm=T))
  bias_lower<-c(bias_lower,quantile((data$postmean),0.025,na.rm=T))
  
  ESSS <- c(ESSS, mean(data$postess,na.rm=T))
  ESSS_75<- c(ESSS_75,quantile((data$postess),0.75,na.rm=T))
  ESSS_25<-c(ESSS_25,quantile((data$postess),0.25,na.rm=T))
  ESSS_upper<- c(ESSS_upper,quantile((data$postess),0.975,na.rm=T))
  ESSS_lower<-c(ESSS_lower,quantile((data$postess),0.025,na.rm=T))
  
}


postvar_diff <- postvar-postvar[1]
case <- c('random averaging','random','1 cluster','3 cluster','baseline: iMEM q=10','baseline: iMEM q=5', 'evenly distributed', 'baseline: kmeans k=10','baseline: kmeans k=5','top 5 + combine 5', 'combine top k')
order <- c('baseline: kmeans k=5','baseline: kmeans k=10','baseline: iMEM q=5','baseline: iMEM q=10','combine top k','evenly distributed','top 5 + combine 5', 'random', 'random averaging','3 cluster','1 cluster')

x <- c()
for (k in 1:length(case)){
  w <- which(order==case[k])
  plugin <- 1+(w-1)*1.25
  x <- c(x, plugin)
}

df2 <- as.data.frame(cbind(postvar_diff,postvar,postvar_upper,postvar_lower,postvar_75,postvar_25,RMSE,RMSE_upper,RMSE_lower,RMSE_75,RMSE_25,bias,bias_upper,bias_lower,bias_75,bias_25,ESSS,ESSS_upper,ESSS_lower,ESSS_75,ESSS_25,x))


s2var <- ggplot(data=df2,aes(x=x))+
  geom_boxplot(aes(ymin=postvar_lower,ymax=postvar_upper,lower=postvar_25,upper=postvar_75,middle=postvar,group=x),
               stat = 'identity', width = 0.8)+
  scale_x_continuous(breaks=seq(1,13.5,1.25),labels=c('kmeans k=5','kmeans k=10','iMEM q=5','iMEM q=10','ordered combine', 'evenly distributed', 'top half & combine half', 'random','random averaging','3 cluster','1 cluster')) +
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x))) +
  ylab('Posterior var')+xlab('(b) Scenario II')+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle=30,hjust=1))

s2bias <- ggplot(data=df2,aes(x=x))+
  geom_boxplot(aes(ymin=bias_lower,ymax=bias_upper,lower=bias_25,upper=bias_75,middle=bias,group=x),
               stat = 'identity', width = 0.8)+
  scale_x_continuous(breaks=seq(1,13.5,1.25),labels=c('kmeans k=5','kmeans k=10','iMEM q=5','iMEM q=10','ordered combine', 'evenly distributed', 'top half & combine half', 'random','random averaging','3 cluster','1 cluster')) +ylab('Bias')+xlab('(b) Scenario II')+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle=30,hjust=1))

s2RMSE <- ggplot(data=df2,aes(x=x))+
  geom_boxplot(aes(ymin=RMSE_lower,ymax=RMSE_upper,lower=RMSE_25,upper=RMSE_75,middle=RMSE,group=x),
               stat = 'identity', width = 0.8)+
  scale_x_continuous(breaks=seq(1,13.5,1.25),labels=c('kmeans k=5','kmeans k=10','iMEM q=5','iMEM q=10','ordered combine', 'evenly distributed', 'top half & combine half', 'random','random averaging','3 cluster','1 cluster')) +
  # scale_y_continuous(trans = 'log10',
  #                       breaks = trans_breaks('log10', function(x) 10^x),
  #                       labels = trans_format('log10', math_format(10^.x))) +
  ylab('RMSE')+xlab('(b) Scenario II')+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle=30,hjust=1))

s2esss <- ggplot(data=df2,aes(x=x))+
  geom_boxplot(aes(ymin=ESSS_lower,ymax=ESSS_upper,lower=ESSS_25,upper=ESSS_75,middle=ESSS,group=x),
               stat = 'identity', width = 0.8)+
  scale_x_continuous(breaks=seq(1,13.5,1.25),labels=c('kmeans k=5','kmeans k=10','iMEM q=5','iMEM q=10','ordered combine', 'evenly distributed', 'top half & combine half', 'random','random averaging','3 cluster','1 cluster')) + ylab('ESSS')+xlab('(b) Scenario II')+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle=30,hjust=1))

s2var
s2bias
s2RMSE
s2esss

################################# scenario III #######################################
dir.path <- "../data/Simulation results/Scenario III/"
library(plyr)
filename <- list.files(path = dir.path,pattern = '.csv')
detach(package:plyr)

postvar <- c()
postvar_upper<-c()
postvar_lower<-c()
postvar_75<-c()
postvar_25<-c()
RMSE <-c()
bias <- c()
RMSE_75 <- c()
RMSE_25 <- c()
RMSE_upper <- c()
RMSE_lower <- c()
bias_75 <- c()
bias_25 <- c()
bias_upper <- c()
bias_lower <- c()
ESSS <- c()
ESSS_75 <- c()
ESSS_25 <- c()
ESSS_upper <- c()
ESSS_lower <- c()


for(i in 1:length(filename)) {
  data <- read.csv(paste0(dir.path,filename)[i])
  
  postvar <- c(postvar,median(data$postvar,na.rm=T))
  postvar_upper<-c(postvar_upper,quantile(data$postvar,0.975,na.rm=T))
  postvar_lower<-c(postvar_lower,quantile(data$postvar,0.025,na.rm=T))
  postvar_75<-c(postvar_75,quantile(data$postvar,0.75,na.rm=T))
  postvar_25<-c(postvar_25,quantile(data$postvar,0.25,na.rm=T))
  
  RMSE_upper<- c(RMSE_upper,quantile(sqrt((data$postmean)^2+data$postvar),0.975,na.rm=T))
  RMSE_lower<-c(RMSE_lower,quantile(sqrt((data$postmean)^2+data$postvar),0.025,na.rm=T))
  RMSE_75<- c(RMSE_75,quantile(sqrt((data$postmean)^2+data$postvar),0.75,na.rm=T))
  RMSE_25<-c(RMSE_25,quantile(sqrt((data$postmean)^2+data$postvar),0.25,na.rm=T))
  RMSE <- c(RMSE,median(sqrt((data$postmean)^2+data$postvar),na.rm=T))
  
  bias <- c(bias,median(data$postmean,na.rm=T))
  bias_75<- c(bias_75,quantile((data$postmean),0.75,na.rm=T))
  bias_25<-c(bias_25,quantile((data$postmean),0.25,na.rm=T))
  bias_upper<- c(bias_upper,quantile((data$postmean),0.975,na.rm=T))
  bias_lower<-c(bias_lower,quantile((data$postmean),0.025,na.rm=T))
  
  ESSS <- c(ESSS, mean(data$postess,na.rm=T))
  ESSS_75<- c(ESSS_75,quantile((data$postess),0.75,na.rm=T))
  ESSS_25<-c(ESSS_25,quantile((data$postess),0.25,na.rm=T))
  ESSS_upper<- c(ESSS_upper,quantile((data$postess),0.975,na.rm=T))
  ESSS_lower<-c(ESSS_lower,quantile((data$postess),0.025,na.rm=T))
  
}


postvar_diff <- postvar-postvar[1]
case <- c('dMEM 0.1', 'iMEM 0.1', 'dMEM 0.2', 'iMEM 0.2', 'dMEM 0.3', 'iMEM 0.3', 'dMEM 0.4', 'iMEM 0.4', 'dMEM 0.5', 'iMEM 0.5')
method <- rep(c('dMEM', 'iMEM'), 5)
test.prop <- rep(c(1:5)/10, each=2)

order <- case

x <- c()
for (k in 1:length(case)){
  w <- which(order==case[k])
  if(method[k] == 'iMEM'){
    plugin <- 1+(w-1)*1.25
  } else{
    plugin <- 1.25+(w-1)*1.25
  }
  x <- c(x, plugin)
}

df4 <- as.data.frame(cbind(postvar_diff,postvar,postvar_upper,postvar_lower,postvar_75,postvar_25,RMSE,RMSE_upper,RMSE_lower,RMSE_75,RMSE_25,bias,bias_upper,bias_lower,bias_75,bias_25,ESSS,ESSS_upper,ESSS_lower,ESSS_75,ESSS_25,x,test.prop))

s4var <- ggplot(data=df4,aes(x=x,fill=method))+
  geom_boxplot(aes(ymin=postvar_lower,ymax=postvar_upper,lower=postvar_25,upper=postvar_75,middle=postvar,group=x),
               stat = 'identity', width = 0.8)+
  geom_vline(xintercept=seq(3,11,2.5),linetype=2)+
  scale_x_continuous(breaks=seq(1.75,12,2.5),labels=c('10%', '20%','30%','40%','50%')) +
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x))) +
  scale_fill_grey(start = 0.7,end=1) +
  ylab('Posterior var')+xlab('\n (c) Scenario III')+ labs(fill = "Method") +
  theme(legend.position = c(0.92, 0.87),
        legend.background = element_rect(fill = "white", color = "black"))

s4bias <- ggplot(data=df4,aes(x=x,fill=method))+
  geom_boxplot(aes(ymin=bias_lower,ymax=bias_upper,lower=bias_25,upper=bias_75,middle=bias,group=x),
               stat = 'identity', width = 0.8)+
  geom_vline(xintercept=seq(3,11,2.5),linetype=2)+
  scale_x_continuous(breaks=seq(1.75,12,2.5),labels=c('10%', '20%','30%','40%','50%')) +
  scale_fill_grey(start = 0.7,end=1) +
  ylab('Bias')+xlab('\n (c) Scenario III')+labs(fill = "Method") +
  theme(legend.position = c(0.92, 0.87),
        legend.background = element_rect(fill = "white", color = "black"))

s4RMSE <- ggplot(data=df4,aes(x=x,fill=method))+
  geom_boxplot(aes(ymin=RMSE_lower,ymax=RMSE_upper,lower=RMSE_25,upper=RMSE_75,middle=RMSE,group=x),
               stat = 'identity', width = 0.8)+
  geom_vline(xintercept=seq(3,11,2.5),linetype=2)+
  scale_x_continuous(breaks=seq(1.75,12,2.5),labels=c('10%', '20%','30%','40%','50%')) +
  scale_fill_grey(start = 0.7,end=1) +
  ylab('RMSE')+xlab('\n (c) Scenario III')+labs(fill = "Method") +
  theme(legend.position = c(0.92, 0.87),
        legend.background = element_rect(fill = "white", color = "black"))

s4esss <- ggplot(data=df4,aes(x=x,fill=method))+
  geom_boxplot(aes(ymin=ESSS_lower,ymax=ESSS_upper,lower=ESSS_25,upper=ESSS_75,middle=ESSS,group=x),
               stat = 'identity', width = 0.8)+
  geom_vline(xintercept=seq(3,11,2.5),linetype=2)+
  scale_x_continuous(breaks=seq(1.75,12,2.5),labels=c('10%', '20%','30%','40%','50%')) + 
  scale_fill_grey(start = 0.7,end=1) +
  ylab('ESSS')+xlab('\n (c) Scenario III')+labs(fill = "Method") +
  theme(legend.position = 'bottom',
        legend.background = element_rect(fill = "white", color = "black"))

s4var
s4bias
s4RMSE
s4esss

################################# Combine #######################################
s4legend <-  ggplot(data=df4,aes(x=x,fill=method))+
  geom_boxplot(aes(ymin=postvar_lower,ymax=postvar_upper,lower=postvar_25,upper=postvar_75,middle=postvar,group=x),
               stat = 'identity', width = 0.8)+
  scale_fill_grey(start = 0.7,end=1) +
  labs(fill = "Method") +
  theme(legend.position = 'bottom',
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.key.height= unit(0.15, 'inch'),
        legend.key.width= unit(0.35, 'inch'),
        legend.background = element_rect(fill = "white", color = "black"))

temp <- ggplotGrob(s4legend)
leg_index <- which(sapply(temp$grobs, function(x) x$name) == "guide-box")
legend <- temp$grobs[[leg_index]]

s1var.plot <- s1var + theme(legend.position="none",axis.title.x = element_blank())
s2var.plot <- s2var + theme(legend.position="none",axis.title.x = element_blank())
s4var.plot <- s4var + theme(legend.position="none",axis.title.x = element_blank())


s1bias.plot <- s1bias + theme(legend.position="none",axis.title.x = element_blank())
s2bias.plot <- s2bias + theme(legend.position="none",axis.title.x = element_blank())
s4bias.plot <- s4bias + theme(legend.position="none",axis.title.x = element_blank())


s1RMSE.plot <- s1RMSE + theme(legend.position="none",axis.title.x = element_blank())
s2RMSE.plot <- s2RMSE + theme(legend.position="none",axis.title.x = element_blank())
s4RMSE.plot <- s4RMSE + theme(legend.position="none",axis.title.x = element_blank())

s1esss.plot <- s1esss + theme(legend.position="none")
s2esss.plot <- s2esss + theme(legend.position="none")
s4esss.plot <- s4esss + theme(legend.position="none",axis.title.x = element_blank())


grid.arrange(arrangeGrob(s1var.plot,s1bias.plot,s1RMSE.plot,
                         nrow=1,heights=2.9,
                         bottom = '(a) Scenario I
                         '),
             arrangeGrob(s2var.plot,s2bias.plot,s2RMSE.plot,
                         nrow=1,heights=2.9,
                         bottom = '(b) Scenario II
                         '),
             arrangeGrob(s4var.plot,s4bias.plot,s4RMSE.plot,legend,
                         layout_matrix=rbind(c(1,2,3),c(NA,4,NA)),
                         nrow=2,heights=c(2.9,0.4),
                         bottom = '\n (c) Scenario III
                         '),
             heights=c(2.9,2.9,3.2),
             nrow=3)

grid.arrange(arrangeGrob(s1esss.plot,s2esss.plot,
                         nrow=1,heights=2.9),
             arrangeGrob(s4esss.plot,legend,
                         layout_matrix=rbind(c(NA,1,NA),c(NA,2,NA)),
                         nrow=2,heights=c(2.9,0.4),widths=c(0.8,2.5,0.8)),
             heights=c(2.9,3.2),
             nrow=2,
             bottom = '(c) Scenario III
                         ')
