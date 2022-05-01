
library(raster)
library(dismo)
library(gbm)
#library(mgcv)


## cross-validation individual models, only for MW2-R25
## compare climate only to climate + biotic
load("InputData/Climate5kmfootprintv2.Rdata")

sp.list <- gsub('.siteIntersAllWebs.Rdata','',list.files(path = "IndivModels",pattern="siteIntersAll")) #367 

sp.list <- sp.list[341:length(sp.list)]

source('Functions/FunctionsBRT.R')



#run models
#runIndModBRTCrossVal(sp.list=sp.list,clim=clim,bio1="MW2Prey",bio2="MW2CompR25",type="R25",metawebName="MW2")


# run models with spatial north-south block and west-east block
runIndModBRTCrossValSpatial(sp.list=sp.list,clim=clim,bio1="MW2Prey",bio2="MW2CompR25",type="R25",metawebName="MW2")








## calculate predictive performance for model fit to test data (50% split)_____________

library(ROCR)
source('Functions/FunctionsCrossVali.R')


#aggregate all three CV models (all, biotic only, env only) for each species to calculate performance metrics

##bring in validation results
pth <- "BRTs"
prefx <- c('SingleSp6predsCVtest.MW2','SingleSp6predsCVtestBioticOnly.MW2','SingleSp6predsCVtestEnvOnly.MW2',
           'SingleSp6predsCVtestSpatialNorth.MW2','SingleSp6predsCVtestBioticOnlySpatialNorth.MW2','SingleSp6predsCVtestEnvOnlySpatialNorth.MW2',
           'SingleSp6predsCVtestSpatialEast.MW2','SingleSp6predsCVtestBioticOnlySpatialEast.MW2','SingleSp6predsCVtestEnvOnlySpatialEast.MW2')


filz <- list.files(path=pth,pattern=prefx[1])  
sp <- gsub(paste(prefx[1],'.R25.SpSelect.',sep=""),'',filz)
sp <- gsub('.Rdata','',sp)
   


## run function to summarize all test validation data
AUC <- rocobBRT(prefx,metric="auc",sp=sp)

PRC <- rocobBRT(prefx,metric="aucpr",sp=sp)

nms <- gsub("SingleSp6predsCVtes",'',prefx)
colnames(auc) <- nms
colnames(prc) <- nms


saveRDS(AUC, file="AUCperspeciesBRT.Rdata")
saveRDS(PRC, file="AUCPRperspeciesBRT.Rdata")




#-------------------------------------------------------------
#plot performance..
pdf(file="PerformanceMeans95CI_BRTs.pdf",height=5,width=9)


sauc <- na.omit(readRDS("AUCperspeciesBRT.Rdata"))
sprc <- na.omit(readRDS("AUCPRperspeciesBRT.Rdata"))

#plot mean and 95% CI for model performance
mets <- list(sauc,sprc)
labs <- c("AUC","AUPRC")
par(mfrow=c(1,2))
for (i in 1:2) {
  m <- mets[[i]]
  a <- m
  #order by model type (full, biotic only, env only)
  a <- a[,c(1,7,4,2,8,5,3,9,6)]
  
  sds <- apply(a,2,sd)
  means <- apply(a,2,mean)
  
  er <- qnorm(0.975)*(sds/sqrt(nrow(a)))
  lw <- means - er
  hg <- means + er
  
  
  plot(c(means[1:3],NA,means[4:6],NA,means[7:9]),ylim=c(0.6,1),pch=19,col=c('azure3','darkkhaki','darkolivegreen',NA),xaxt='n',
       ylab=labs[i],xlab="",cex=1.8,cex.axis=1.2)
  arrows(1:11,c(lw[1:3],NA,lw[4:6],NA,lw[7:9]),1:11,c(hg[1:3],NA,hg[4:6],NA,hg[7:9]),
         code=3,length=0,angle=90,lwd=2,col=c('azure3','darkkhaki','darkolivegreen',NA))
  axis(1,at=c(2,6,10),labels=c('Full Model','Biotic Only','Env. Only'),cex.axis=1.2)
  if (i==2)
    legend("topright",       
           legend = c("Random", "Block - Longitude","Block - Latitude"),
           col = c('azure3','darkkhaki','darkolivegreen'),
           pch = 19,cex=1.2)
}
dev.off()








