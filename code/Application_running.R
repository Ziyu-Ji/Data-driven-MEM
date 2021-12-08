
############################### Packages, data and functions ############################
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtable)
library(gridExtra)
library(ggtext)
library(grid)

# data and function loaded
source('Application_functions.R')
load('../data/realdata.rdata')
realdata <- imem_realdata
rm(imem_realdata)
primid <- unique(realdata$surveycode)

###################### Application results for different dMEM clustering methods#########################
## topk works the best
gendat.emotbymode(primid,mtd='even',data=realdata)
gendat.emotbymode(primid,mtd='topk',data=realdata)
gendat.emotbymode(primid,mtd='random',data=realdata)
gendat.emotbymode(primid,mtd='random',m=10,data=realdata)

# QC check table
DMEM_qc(primid,data=realdata)

#################################### Performance analysis########################################
qc <- read.csv('../data/Application results/app_qc_all.csv',header=T)
topk <- read.csv('../data/Application results/application_1m_topkmethod.csv',header=T)
topk$mode <- substr(topk$mode,1,nchar(topk$mode)-7)

## remove 'IN_VEHICLE'
topk <- filter(topk, mode != 'IN_VEHICLE')
qc <- filter(qc, mode != 'IN_VEHICLE')

## general amount of difference between iMEM and dMEM
topk$diff <- topk$dmem.pctdecsd - topk$imem.pctdecsd
t.test(topk$diff)

## ESSS of iMEM and dMEM
mean(topk$imem.esss)
mean(topk$dmem.esss)

## mean num of supple for car = 75.50
mean(qc[which(qc$mode=='CAR'),1])
## mean num of supple for walk = 48.10
mean(qc[which(qc$mode=='WALK'),1])
## mean num of supple for bus = 6.20
mean(qc[which(qc$mode=='BUS'),1])
## mean num of supple for bike = 6.61
mean(qc[which(qc$mode=='BIKE'),1])


# hist plot for percentage of selected sources
hist(qc[,1]/qc[,2])
mean(qc[which(qc$mode=='CAR'),1]/qc[which(qc$mode=='CAR'),2])
mean(qc[which(qc$mode=='WALK'),1]/qc[which(qc$mode=='WALK'),2])
mean(qc[which(qc$mode=='BUS'),1]/qc[which(qc$mode=='BUS'),2])
mean(qc[which(qc$mode=='BIKE'),1]/qc[which(qc$mode=='BIKE'),2])


# when imem performs bad (improve < 20% => top 93)
## how many of them have dmem performance better than imem? 87/93=93.55%
length(which(topk$imem.pctdecsd<0.2))
dmem_perf <- topk$dmem.pctdecsd[order(topk$imem.pctdecsd,decreasing=F)]
dmem.imem.diff <- dmem_perf[1:93] - sort(topk$imem.pctdecsd,decreasing=F)[1:93]
length(which(dmem.imem.diff>0))
plot(x=sort(topk$imem.pctdecsd,decreasing=F)[1:93],y=dmem.imem.diff)

## how much improvement? mean = 0.24
hist(dmem.imem.diff,breaks=20)
t.test(dmem.imem.diff)

# when dmem performs bad (improve < 20% => top 68), how many of them have imem performance better than dmem? 7/68=10.29%
length(which(topk$dmem.pctdecsd<0.2))
imem_perf <- topk$imem.pctdecsd[order(topk$dmem.pctdecsd,decreasing=F)]
imem.dmem.diff <- imem_perf[1:68] - sort(topk$dmem.pctdecsd,decreasing=F)[1:68]
length(which(imem.dmem.diff>0))
plot(x=sort(topk$dmem.pctdecsd,decreasing=F)[1:68],y=imem.dmem.diff)

## how much improvement? mean = -0.136
hist(imem.dmem.diff)
t.test(imem.dmem.diff)

## merge the dataset for plotting
qc$index <- paste0(qc$mode,qc$outcm,qc$primid)
topk$index <- paste0(topk$mode,topk$outcm,topk$primid)
qc.topk <- merge(qc, topk, by.x='index', by.y='index')

###################################### Plots #############################################

# histograms on results
hist(qc.topk$imem.pctdecsd,breaks = 30)
hist(qc.topk$dmem.pctdecsd,breaks = 30)

hist(qc.topk$dmem.esss,breaks = 30)
hist(qc.topk$imem.esss,breaks = 30)

imprv <- ggplot(data=qc.topk)+ stat_density(aes(x=imem.pctdecsd*100, y=..density..,linetype="dashed"), geom="line",position="identity")+ stat_density(aes(x=dmem.pctdecsd*100,y=..density..,linetype="solid"), geom="line",position="identity") + xlim(-50,100) +
  xlab('SD reduction (%)') + ylab('Density') +
  scale_linetype_manual(name="Method", values=c(dashed="dashed", solid="solid"), labels=c('iMEM', 'dMEM')) 

esss <- ggplot(data=qc.topk)+ stat_density(aes(x=imem.esss, y=..density..,linetype="dashed"), geom="line",position="identity")+ stat_density(aes(x=dmem.esss,y=..density..,linetype="solid"), geom="line",position="identity") + xlim(0,5000) +
  xlab('ESSS') + ylab('Density') +
  scale_linetype_manual(name="Method", values=c(dashed="dashed", solid="solid"), labels=c('iMEM', 'dMEM')) 

diff_hist <- ggplot(data=qc.topk,aes(x=diff*100, y=..density..))+ geom_histogram(bins=100) + geom_density()+ xlim(-50,50) +
  xlab('Pairwise difference in SD reduction (%)') + ylab('Density')

grid.arrange(imprv,diff_hist,esss,widths=5,heights=c(3,3,3),nrow=3)


# number of selected sources
num_source_data <- qc.topk %>%
  group_by(mode.x) %>%
  summarise(mean_select = round(mean(selected.supple),digits=1),
            source_95 = round(quantile(selected.supple,0.95),digits=1),
            source_5 = round(quantile(selected.supple,0.05),digits=1),
            mean_total = round(mean(total.supple),digits=1))

num.select.plot.label <- num_source_data$mode.x
num_source_data$mode.x <- c(0,3.5,7,10.5)

num.select.plot <- 
  ggplot(data=num_source_data) +
  geom_bar(aes(x=mode.x, y=mean_total, fill = "grey30"), stat = 'identity', width = 0.8) +
  geom_bar(aes(x=mode.x+1, y=mean_select, fill = "grey50"), stat = 'identity', width = 0.8) +
  geom_bar(aes(x=mode.x+2, y=10, fill = "grey70"), stat = 'identity', width = 0.8) +
  geom_errorbar(aes(x=mode.x+1, ymin=source_5, ymax=source_95), width=0.2) + ylab('Number of selected supplementary sources') +
  scale_y_continuous(breaks = seq(0, 260, by = 20)) +
  scale_x_continuous(breaks = c(1,4.5,8,11.5),label=num.select.plot.label) +
  scale_fill_manual(values = c('grey70','grey50','grey30'), breaks = c('grey70','grey50','grey30'), label = c('iMEM','dMEM','Total'), name = 'Method') +
  theme(axis.title.y=element_blank(),axis.text = element_text(size=10), axis.title.x = element_text(size=12)) +
  coord_flip() 


# examples of marginal weight plots
example <- qc.topk[c(21,97,307,1828),]

exp1 <- ggplot() + geom_point(aes(x=c(1:length(as.numeric(unlist(strsplit(example$weights[1], ","))))),
                                  y=as.numeric(unlist(strsplit(example$weights[1], ",")))), shape = 1, size = 3) + 
  geom_vline(aes(xintercept =example$selected.supple[1], linetype="solid"), show.legend=TRUE) + 
  geom_vline(aes(xintercept = 10, linetype="dashed"), show.legend=TRUE) + 
  scale_linetype_manual(name="Selected sources", values=c(dashed="dashed", solid="solid"), labels=c('iMEM', 'dMEM')) +
  labs(title='ID 8030, Trip mode: BIKE, Emotion: meaningful', caption = '(SD reduction: iMEM 37.6%, dMEM 76.0%)') + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(size=12),plot.caption = element_text(size=10),axis.text = element_text(size=10))

exp2 <- ggplot() + geom_point(aes(x=c(1:length(as.numeric(unlist(strsplit(example$weights[2], ","))))),
                                  y=as.numeric(unlist(strsplit(example$weights[2], ",")))), shape = 1, size = 3) + 
  geom_vline(xintercept =example$selected.supple[2]) + 
  geom_vline(xintercept = 10, linetype="dashed") + 
  labs(title='ID 4097, Trip mode: BUS, Emotion: pain', caption='(SD reduction: iMEM -65.6%, dMEM 32.1%)') + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(size=12),plot.caption = element_text(size=10),axis.text = element_text(size=10))

exp3 <- ggplot() + geom_point(aes(x=c(1:length(as.numeric(unlist(strsplit(example$weights[3], ","))))),
                                  y=as.numeric(unlist(strsplit(example$weights[3], ",")))), shape = 1, size = 3) + 
  geom_vline(xintercept =example$selected.supple[3]) + 
  geom_vline(xintercept = 10, linetype="dashed") + 
  scale_x_continuous(breaks = seq(0, 300, by = 30)) +
  labs(title='ID 5088, Trip mode: CAR, Emotion: happy', caption = '(SD reduction: iMEM 58.1%, dMEM 86.8%)') + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(size=12),plot.caption = element_text(size=10),axis.text = element_text(size=10))

exp4 <- ggplot() + geom_point(aes(x=c(1:length(as.numeric(unlist(strsplit(example$weights[4], ","))))),
                                  y=as.numeric(unlist(strsplit(example$weights[4], ",")))), shape = 1, size = 3) + 
  geom_vline(xintercept =example$selected.supple[4]) + 
  geom_vline(xintercept = 10, linetype="dashed") + 
  scale_x_continuous(breaks = seq(0, 150, by = 20)) +
  labs(title='ID 5056, Trip mode: WALK, Emotion: tired', caption= '(SD reduction: iMEM 39.3%, dMEM 72.4%)') + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(size=12),plot.caption = element_text(size=10),axis.text = element_text(size=10))

legend <- gtable_filter(ggplotGrob(exp1), "guide-box")

exp1 <- exp1 + theme(legend.position="none")
exp2 <- exp2 + theme(legend.position="none")
exp3 <- exp3 + theme(legend.position="none")
exp4 <- exp4 + theme(legend.position="none")

label <- textGrob("Marginal weight", rot = 90, vjust = 0.5)


weight.dist.plot <- grid.arrange(label,
                                 arrangeGrob(exp1, exp2, exp3, exp4,
                                             nrow=2,
                                             bottom = 'Sorted supplementary sources'),
                                 legend,nrow=1,
                                 widths=c(0.2,3.5,0.6), 
                                 heights=6,
                                 bottom='(b) Examples of sorted marginal weights and selected sources
  ')

legend1 <- gtable_filter(ggplotGrob(num.select.plot), "guide-box")

num.select.plot <- num.select.plot + theme(legend.position="none")

label1 <- textGrob("Trip mode", rot = 90, vjust = 0.5)

num.select.plot.fin <- grid.arrange(label1,num.select.plot,legend1,
                                    heights=0.5, nrow=1, widths=c(0.1,3.5,0.6), 
                                    bottom='(a) Bar blot for the numbers of selected supplementary sources
  
  ')

grid.arrange(num.select.plot.fin,weight.dist.plot,nrow=2,heights=c(3.5,5),layout_matrix=rbind(c(1,1,1,1,1),c(2,2,2,2,2)))

