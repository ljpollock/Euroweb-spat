


## cross-validation individual models, only for MW2-R25
## compare climate only to climate + biotic

reticulate::use_condaenv("r-tensorflow", required = TRUE)


##### Load packages and data ####
library(greta)
library(tidyr)
library(dplyr)
library(bayesplot)
library(greta)
library(igraph)
library(DiagrammeR)

library(Matrix)
library(reshape2)


source("Functions/FunctionsPythonModels.R")
source('FunctionsCrossVali.R')

#climate inputs
load("InputData/Climate5kmfootprintv2.Rdata")

# dir for individual mods
#Analyses3Feb20/IndivModels")

sp.list <- gsub('.siteIntersAllWebs.Rdata','',list.files(pattern="siteIntersAll")) #367 



#run models - Careful- runs ALL models
#runIndModPythonCrossVal(sp.list=sp.list,clim=clim,bio1="MW2Prey",bio2="MW2CompR25",type="R25",metawebName="MW2")


# run models with spatial north-south block and west-east block
#runIndModPythonCrossValSpatial(sp.list=sp.list,clim=clim,bio1="MW2Prey",bio2="MW2CompR25",type="R25",metawebName="MW2")



## calculate predictive performance for model fit to test data (50% split)_____________

library(ROCR)

#aggregate all three CV models (all, biotic only, env only) for each species to calculate performance metrics


##bring in validation results
prefx <- c('SingleSp6predsCVtest.MW2','SingleSp6predsCVtestBioticOnly.MW2','SingleSp6predsCVtestEnvOnly.MW2',
           'SingleSp6predsCVtestSpatialNorth.MW2','SingleSp6predsCVtestBioticOnlySpatialNorth.MW2','SingleSp6predsCVtestEnvOnlySpatialNorth.MW2',
           'SingleSp6predsCVtestSpatialEast.MW2','SingleSp6predsCVtestBioticOnlySpatialEast.MW2','SingleSp6predsCVtestEnvOnlySpatialEast.MW2')

val <- lapply(prefx,agg,n=length(list.files(pattern=prefx)))
names(val) <- gsub("SingleSp6predsCVtes",'',prefx)


indx <- Reduce(intersect, list(rownames(val$t.MW2),rownames(val$tSpatialEast.MW2),
                               rownames(val$tSpatialNorth.MW2)))

valint <- lapply(seq_along(val), function(i) {
  dat <- val[[i]]
  subset(dat, row.names(dat) %in% indx)
})
names(valint) <- gsub("SingleSp6predsCVtes",'',prefx)

#select probabilities from list
valpr <- lapply(valint, '[[', 3)
#select occurrence 
valoccur <- lapply(valint, '[[', 2)


#split species sites
l <- data.frame(do.call('rbind', strsplit(as.character(rownames(valint[[1]])),'.',fixed=TRUE)))
colnames(l) <- c("Species","Site")

sp.long <- l$Species

#need species in long form, occurrences and probs from validation plots
auc <- rocob(splong=sp.long,occur=valoccur,pr=valpr,metric='auc')
prc <- rocob(splong=sp.long,occur=valoccur,pr=valpr,metric='aucpr')

colnames(auc) <- names(valint)
colnames(prc) <- names(valint)


save(auc, file="ResultFiles/AUCperspecies.Rdata")
save(prc, file="ResultFiles/AUCPRperspecies.Rdata")


#---------------------------------------------------------------------------------------



# calculate species median latitude.. 
load('InputData/SiteBySpecies21Feb.Rdata')
load("InputData/Climate5kmfootprintv2.Rdata")


load("InputData/SpeciesLevelWebStatisticsMW2R25.Rdata")
spl <- do.call("cbind",SpWebStats)

lat <- clim[,'y']
names(lat) <- row.names(clim)

sp.list <- row.names(spl[which(spl[,3]>0),])

median.lat <- data.frame(matrix(NA,ncol=1,nrow=length(sp.list)))
rownames(median.lat) <- sp.list
colnames(median.lat) <- "medianLat"

for(i in 1:length(sp.list)) {
  sp <- sp.list[i]
  sites <- site.by.species[,sp]
  present <- sites[which(sites==1)]
  m.lat <- median(  lat[match(names(present),names(lat))] ,na.rm=T)
  median.lat[sp,] <- m.lat
}

save(median.lat, file="ResultFiles/speciesMedianLatitude.Rdata")

#----------------------------------------------------------------------------------------





#---------------------------------------------------------------------------------------
#  combine with other species-level metrics and plot results
# bring in species taxonomy
taxo <- read.csv(file="InputData/TetrapodsIDSpecies.csv")

#median latitude
load("ResultFiles/speciesMedianLatitude.Rdata")
#auc and prc stats
load("ResultFiles/AUCperspecies.Rdata")
load("ResultFiles/AUCPRperspecies.Rdata")


#network metrics summarized for all species
load("InputData/SpeciesLevelWebStatisticsMW2R25.Rdata")

# with network data from Braga 2019
load("InputData/DF.Rdata")



## merge datasets into a single species-level 

spl <- do.call("cbind",SpWebStats)

spdat <- merge(taxo,DF[,c('TL','Eig.cent')],by.x='ID',by.y=0,all.x=T,all.y=F)

spdat <- merge(spdat,median.lat,by.x='ID',by.y=0,all.x=T,all.y=F)

spdat <- merge(spdat,spl,by.x='ID',by.y=0,all=F)


sauc <- merge(spdat,auc,by.x='ID',by.y=0,all=F)
sprc <- merge(spdat,prc,by.x='ID',by.y=0,all=F)
sauc <- na.omit(sauc)
sprc <- na.omit(sprc)



#-------------------------------------------------------------
#plot performance..
pdf(file="PerformanceMeans95CI.pdf",height=5,width=9)

#plot mean and 95% CI for model performance
mets <- list(sauc,sprc)
labs <- c("AUC","AUPRC")
par(mfrow=c(1,2))
for (i in 1:2) {
m <- mets[[i]]
a <- m[,18:ncol(m)]
#order by model type (full, biotic only, env only)
a <- a[,c(1,7,4,2,8,5,3,9,6)]

sds <- apply(a,2,sd)
means <- apply(a,2,mean)

er <- qnorm(0.975)*(sds/sqrt(nrow(a)))
lw <- means - er
hg <- means + er


plot(c(means[1:3],NA,means[4:6],NA,means[7:9]),ylim=c(0.45,1),pch=19,col=c('azure3','darkkhaki','darkolivegreen',NA),xaxt='n',
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



#-===============================================================
#  Scatterplot of AUC for all species for appendix

BirdCol <- rgb(col2rgb("#7570b3")[1],col2rgb("#7570b3")[2],col2rgb("#7570b3")[3],max = 255, alpha = 90)
MamCol <- rgb(col2rgb("#d95f02")[1],col2rgb("#d95f02")[2],col2rgb("#d95f02")[3],max = 255, alpha = 90)
AmphCol <- rgb(col2rgb("#e7298a")[1],col2rgb("#e7298a")[2],col2rgb("#e7298a")[3],max = 255, alpha = 90)
ReptCol <- rgb(col2rgb("#1b9e77")[1],col2rgb("#1b9e77")[2],col2rgb("#1b9e77")[3],max = 255, alpha = 90)



pdf(file="PerformanceScatterPlotBioVEnv.pdf",height=7,width=9)

mets <- list(sauc,sprc)
labs <- c("AUC","AUPRC")
main <- c("Random","")
main2 <- c("Block - Longitude","")
main3 <- c("Block - Latitude","")
par(mfrow=c(2,3))
for (i in 1:2) {
  m <- mets[[i]]
  a <- m[,c(2,18:ncol(m))]
  colz = c(AmphCol, BirdCol, MamCol, ReptCol)[as.factor(a[,1])]
  #order by model type (full, biotic only, env only)
  plot(a$tEnvOnly.MW2,a$tBioticOnly.MW2,xlim=c(-0.01,1.0),ylim=(c(0,1.0)),ylab=paste(labs[i],"Biotic-Only Model"),
       xlab=paste(labs[i],"Env-Only Model"),col=colz,pch=19,main=main[i])
  abline(abline(a=0, b=1,col='grey60'))
  plot(a$tEnvOnlySpatialEast.MW2,a$tBioticOnlySpatialEast.MW2,xlim=c(0,1.0),ylim=(c(-0.01,1.0)),ylab=paste(labs[i],"Biotic-Only Model"),
       xlab=paste(labs[i],"Env-Only Model"),col=colz,pch=19,main=main2[i])
  abline(abline(a=0, b=1,col='grey60'))
  plot(a$tEnvOnlySpatialNorth.MW2,a$tBioticOnlySpatialNorth.MW2,xlim=c(0,1.0),ylim=(c(-0.01,1.0)),ylab=paste(labs[i],"Biotic-Only Model"),
       xlab=paste(labs[i],"Env-Only Model"),col=colz,pch=19,main=main3[i])
  abline(abline(a=0, b=1,col='grey60'))
}
dev.off()



#---------------------------------------------------------------------
# boxplots of orders

pdf(file='checknames.pdf',width=10)
boxplot(a[,2]~a[,1],yaxt='n',las=2)
dev.off()



#  Boxplot of AUC for all species orders for appendix
pdf(file="BoxplotsPerfOrdersAUC.pdf",width=15,height=8)
m <- sauc
a <- m[,c(3,18:ncol(m))]
ords  <- m[,c(2:3)][!duplicated(m[,3]),]
l <- factor(a[,1], levels = rev(as.factor(ords$Order)))
colz = c(AmphCol, BirdCol, MamCol, ReptCol)[rev(as.factor(ords$Class))]

par(mfrow=c(1,7))

par(mar=c(7,10,5,1))
boxplot(a[,2]~l,horizontal = TRUE,col=colz,yaxt='n',xlab='',ylab="")
axis(2,at=1:nrow(ords),labels=rev(ords[,2]),las=2)

labs <- rep(c('Full',"Biotic","Env."),2)
labs2 <- c(rep('Random',3),rep('Block - Latitude',3))
x=0
par(mar=c(7,2,5,1))
for (i in 2:7){
  x=x+1
  boxplot(a[,i]~l,horizontal = TRUE,col=colz,yaxt='n',xlab=paste(labs[x],"AUC",labs2[x]),ylab="",ylim=c(0,1))
  abline(v=0.5,col='grey')
}
dev.off()


#  Boxplot of AU PRC for all species orders for appendix
pdf(file="BoxplotsPerfOrdersAUPRC.pdf",width=15,height=8)
m <- sprc
a <- m[,c(3,18:ncol(m))]
ords  <- m[,c(2:3)][!duplicated(m[,3]),]
l <- factor(a[,1], levels = rev(as.factor(ords$Order)))
colz = c(AmphCol, BirdCol, MamCol, ReptCol)[rev(as.factor(ords$Class))]

par(mfrow=c(1,7))

par(mar=c(7,10,5,1))
boxplot(a[,2]~l,horizontal = TRUE,col=colz,yaxt='n',xlab='',ylab="")
axis(2,at=1:nrow(ords),labels=rev(ords[,2]),las=2)

labs <- rep(c('Full',"Biotic","Env."),2)
labs2 <- c(rep('Random',3),rep('Block - Latitude',3))
x=0
par(mar=c(7,2,5,1))
for (i in 2:7){
  x=x+1
  boxplot(a[,i]~l,horizontal = TRUE,col=colz,yaxt='n',xlab=paste(labs[x],"AUPRC",labs2[x]),ylab="",ylim=c(0,1))
  abline(v=0.5,col='grey')
}
dev.off()










#-===============================================================
#  What predicts AUC, AUROC, and the loss of predictive power between bio and clim runs?

#------------------------------------

#AUC
m <- sauc[,c(2,9,10,12,13,18:ncol(sauc))]

params.auc <- vector(mode="list", length=9)
names(params.auc) <- c(colnames(m[,6:ncol(m)]))

for (i in 1:length(names(params.auc))){
    nm <- names(params.auc)[i]
    mod <- glm(m[,nm] ~ log(m$RangeSize) + m$medianLat,family='binomial')
    params.auc[[nm]] <- summary(mod)$coef
}

#AUPRC
rm(m)
m <- sprc[,c(2,9,10,12,13,18:ncol(sprc))]
m[,6:ncol(m)][m[,6:ncol(m)]>1] <- 1

params.prc <- vector(mode="list", length=9)
names(params.prc) <- c(colnames(m[,6:ncol(m)]))

for (i in 1:length(names(params.prc))){
  nm <- names(params.prc)[i]
  mod <- glm(m[,nm] ~ log(m$RangeSize) + m$medianLat,family='binomial')
  params.prc[[nm]] <- summary(mod)$coef
}

# keep only significant effects

a <- do.call("rbind",params.auc)
a <- a[-which(rownames(a)%in%'(Intercept)'),]
rownames(a) <- paste(rep(names(params.auc),each=2),row.names(a),sep='--')
sig.auc <- a[which(a[,'Pr(>|z|)'] < 0.05),]


b <- do.call("rbind",params.prc)
b <- b[-which(rownames(b)%in%'(Intercept)'),]
rownames(b) <- paste(rep(names(params.prc),each=2),row.names(b),sep='--')
sig.prc <- b[which(b[,'Pr(>|z|)'] < 0.05),]

# compare only bio vs env models 
au <- sig.auc[grep('EnvOnly|BioticOnly',rownames(sig.auc)),]

pr <- sig.prc[grep('EnvOnly|BioticOnly',rownames(sig.prc)),]


save(au,file="Code/PerformanceGLM_AUC.R")
save(pr,file="Code/PerformanceGLM_AUPRC.R")

# Fig. plot response of AUC/PRC to Range size ==========================

# cols for random, block-long, block-lat
c1 <- rgb(col2rgb("azure3")[1],col2rgb("azure3")[2],col2rgb("azure3")[3],max = 255, alpha = 90)
c2 <- rgb(col2rgb("darkkhaki")[1],col2rgb("darkkhaki")[2],col2rgb("darkkhaki")[3],max = 255, alpha = 90)
c3 <- rgb(col2rgb("darkolivegreen")[1],col2rgb("darkolivegreen")[2],col2rgb("darkolivegreen")[3],max = 255, alpha = 90)
colz=c(c1,c2,c3)
colzF=c('azure3','darkkhaki','darkolivegreen')

#--
source('Functions/FunctionsCrossVali.R')


pdf(file='Figures/ResponseCurvesPredAUC.pdf',width=8,height=8)
par(mfrow=c(2,2))

#1. AUC 
m <- sauc[,c(2,9,10,12,13,18:ncol(sauc))]

x=log(m[,'RangeSize'])
x2=m[,'medianLat']

pred.mat <- data.frame(cbind( seq(min(x),max(x),length.out=100),rep(mean(x2),100)))
colnames(pred.mat) <- c('x','x2')

#biotic
nmz <- c('tBioticOnly.MW2',
         'tBioticOnlySpatialEast.MW2',
         'tBioticOnlySpatialNorth.MW2')

plot.preds.auc(nmz=nmz,m=m,colz=colz,colzF=colzF,x=x,x2=x2,ylab='AUC (Biotic Only Models)')

#2. AUC env
nmz <- c('tEnvOnly.MW2',
         'tEnvOnlySpatialEast.MW2',
         'tEnvOnlySpatialNorth.MW2')

plot.preds.auc(nmz=nmz,m=m,colz=colz,colzF=colzF,x=x,x2=x2,ylab='AUC (Env. Only Models)')

#3. AUPRC  - biotic----

m <- sprc[,c(2,9,10,12,13,18:ncol(sprc))]
m[,6:ncol(m)][m[,6:ncol(m)]>1] <- 1

x=log(m[,'RangeSize'])
x2=m[,'medianLat']
pred.mat <- data.frame(cbind( seq(min(x),max(x),length.out=100),rep(mean(x2),100)))
colnames(pred.mat) <- c('x','x2')

nmz <- c('tBioticOnly.MW2',
         'tBioticOnlySpatialEast.MW2',
         'tBioticOnlySpatialNorth.MW2')

plot.preds.auc(nmz=nmz,m=m,colz=colz,colzF=colzF,x=x,x2=x2,ylab='AUPRC (Biotic Only Models)')


#4. AUPRC  
nmz <- c('tEnvOnly.MW2',
         'tEnvOnlySpatialEast.MW2',
         'tEnvOnlySpatialNorth.MW2')

plot.preds.auc(nmz=nmz,m=m,colz=colz,colzF=colzF,x=x,x2=x2,ylab='AUPRC (Env. Only Models)')

dev.off()
#----------




#----------------------------------------

# plot betas 

colz = c(AmphCol, BirdCol, MamCol, ReptCol)[as.factor(mod.dat$group)]
t <- c("A","B","M","R") 
colz2 <- c("#e7298a","#7570b3","#d95f02","#1b9e77")


#loss of predictive power as percent of original AUC
biominenv <- (sauc$t.MW2-sauc$tBioticOnly.MW2)/sauc$t.MW2




#-------------------------------------

## mis-calibration of enviro predictors

#-------------------------------------


## 1. gold-standard betas- parameter estimates drawn from full models (MW 2, R25)

setwd("~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/Analyses3Feb20/IndivModels")
files <- list.files(pattern='SingleSp6preds.MW2.R25')

modstat <- lapply(files,function (x) get(load(x)))
names(modstat) <- gsub(".Rdata","",gsub("SingleSp6preds.MW2.R25.SpSelect.","",files))

params.gold <- sapply(modstat, "[",1)
meanparams.gold <- t(as.data.frame(sapply(params.gold, function(x) x[,1])))
rownames(meanparams.gold) <- gsub('.statistics','',rownames(meanparams.gold))
colnames(meanparams.gold) <- c('alpha','prey','comp','b12','b6','b15','foot')

#extract sign and signficance - credible interval entirely less than zero, 
#overlaps zero, or entirely greater than zero for prey and competitor terms
credible <- sapply(modstat, "[",2)
credible.2.5 <- t(as.data.frame(sapply(credible, function(x) x[,1])))
credible.97.5 <- t(as.data.frame(sapply(credible, function(x) x[,5])))
credible.2.5.sign <- sign(credible.2.5)
credible.97.5.sign <- sign(credible.97.5)

cred.sum <- credible.2.5.sign + credible.97.5.sign
cred.sum <- cred.sum[,2:3]
cred.sum[cred.sum==2] <- "Positive"
cred.sum[cred.sum==-2] <- "Negative"
cred.sum[cred.sum==0] <- "Overlaps 0"

cat <- gsub("[0-9]*.quantiles","",rownames(cred.sum))
cred.sum <- cbind(cat,cred.sum) 

dat <- data.frame(cred.sum)
colnames(dat) <- c("Class","Prey","Competitor")

pdf(file="BarplotBetasGoldPercentCredibleIntervals_Prey.pdf")
ggplot(dat) +
  geom_bar(aes(x=Class, fill=Prey,stat="identity"))
dev.off()

pdf(file="BarplotBetasGoldPercentCredibleIntervals_Competitor.pdf")
ggplot(dat) +
  geom_bar(aes(x=Class, fill=Competitor,stat="identity"))
dev.off()



#boxplot of beta MW2R25 for appendix

betas <- merge(spdat,meanparams.gold,by.x='ID',by.y=0,all=F)


#  Boxplot for all species orders for appendix
pdf(file="BoxplotsOrdersBetasGold.pdf",width=15,height=8)
m <- betas
a <- m[,c(3,18:ncol(m))]
ords  <- m[,c(2:3)][!duplicated(m[,3]),]
l <- factor(a[,1], levels = rev(as.factor(ords$Order)))
colz = c(AmphCol, BirdCol, MamCol, ReptCol)[rev(as.factor(ords$Class))]

par(mfrow=c(1,8))

par(mar=c(7,10,5,1))
boxplot(a[,2]~l,horizontal = TRUE,col='white',xaxt='n',yaxt='n',xlab='',ylab="",ylim=c(-4.5,4.5))
axis(2,at=1:nrow(ords),labels=rev(ords[,2]),las=2)

labs <- c("Order","Alpha","Prey","Comp.","Precip.","Min. Temp.","Precip. Season.","Footprint")

par(mar=c(7,2,5,1))
for (i in 2:8){
  boxplot(a[,i]~l,horizontal = TRUE,col=colz,yaxt='n',xlab=labs[i],ylab="",ylim=c(-4.5,4.5))
  abline(v=0.5,col='grey')
}
dev.off()




# compare betas from env-only models to gold standard beta envs

# 2. biotic only - Random, Block lat, block long
meanparams.b <- vector(mode="list", length=3)
names(meanparams.b) <- c('Random','Blong','Blat')

pts <- c('SingleSp6predsCVtrainBioticOnly.MW2.R25.SpSelect.',
         'SingleSp6predsCVtrainBioticOnlySpatialEast.MW2.R25.SpSelect.',
         'SingleSp6predsCVtrainBioticOnlySpatialNorth.MW2.R25.SpSelect.')

for (i in 1:3) {
  s <- list.files(pattern=pts[i])
  comb <- lapply(s,function(x) get(load(x)))
  params <- sapply(comb, "[",1)
  e <- t(as.data.frame(sapply(params, function(x) x[,1])))
  sp <- gsub(pts[i],"",s)
  sp <- gsub(".Rdata","",sp)
  rownames(e) <- sp
  colnames(e) <- c('alpha','prey','comp')
  meanparams.b[[i]] <- e
  rm(s,e,sp)
}

# 3.  env only - Random, Block lat, block long
meanparams.e <- vector(mode="list", length=3)
names(meanparams.e) <- c('Random','Blong','Blat')

pts <- c('SingleSp6predsCVtrainEnvOnly.MW2.R25.SpSelect.',
         'SingleSp6predsCVtrainEnvOnlySpatialEast.MW2.R25.SpSelect.',
         'SingleSp6predsCVtrainEnvOnlySpatialNorth.MW2.R25.SpSelect.')

for (i in 1:3) {
  s <- list.files(pattern=pts[i])
  comb <- lapply(s,function(x) get(load(x)))
  params <- sapply(comb, "[",1)
  e <- t(as.data.frame(sapply(params, function(x) x[,1])))
  sp <- gsub(pts[i],"",s)
  sp <- gsub(".Rdata","",sp)
  rownames(e) <- sp
  colnames(e) <- c('alpha','b12','b6','b15','foot')
  meanparams.e[[i]] <- e
  rm(s,e,sp)
}


## calc differences
## gold - biotic
indx <- Reduce(intersect, list(rownames(meanparams.b$Random),
                               rownames(meanparams.b$Blong),
                               rownames(meanparams.b$Blat),
                               rownames(meanparams.e$Random),
                               rownames(meanparams.e$Blong),
                               rownames(meanparams.e$Blat)))

b <- lapply(seq_along(meanparams.b), function(i) {
  dat <- meanparams.b[[i]]
  subset(dat, row.names(dat) %in% indx)
})

r <- meanparams.gold[rownames(meanparams.gold)%in%indx,c(1:3)]
diff.b.a <- lapply(b, function(x, y) abs(x - y), y = r)
diff.b <- lapply(b, function(x, y) x - y, y = r)


## gold - env
e <- lapply(seq_along(meanparams.e), function(i) {
  dat <- meanparams.e[[i]]
  subset(dat, row.names(dat) %in% indx)
})

r <- meanparams.gold[rownames(meanparams.gold)%in%indx,c(1,4:7)]
diff.e.a <- lapply(e, function(x, y) abs(x - y), y = r)
diff.e <- lapply(e, function(x, y) x - y, y = r)

#ss <- sapply(diff.e, function(x) colSums(x))


#---------------------------------
## plot mean effect (gold standard), bias (mean effect - effect beta CV runs), and miscalibration (abs |mean effect - mean effect|))
pdf(file='BiasMiscalibrationBetas.pdf',width=10,height=5)

par(mar=c(5,7,2,1))

par(mfrow=c(1,3))
a <- meanparams.gold

sds <- apply(a,2,sd)
means <- apply(a,2,mean)

er <- qnorm(0.975)*(sds/sqrt(nrow(a)))
lw <- means - er
hg <- means + er

plot(means[2:7],6:1,pch=19,col='grey35',yaxt='n',cex=1.8,cex.axis=1.2,
     xlim=c(-0.5,1.7),cex.lab=1.2,ylab='',xlab='Mean Effect (Final Model)')

arrows(lw[2:7],6:1,hg[2:7],6:1,
       code=3,length=0,angle=90,lwd=2,col=c('grey35'))
abline(v=0,col='grey')
axis(2,at=6:1,labels=c('Prey','Competitor','Precip.','Min. temp','Prec. seas.','Footprint'),las=1,cex.axis=1.2)


# raw difference (bias)
par(mar=c(5,1,2,1))
for (i in 1:3) {
        b <- diff.b[[i]]
        e <- diff.e[[i]]
        be <- cbind(b[,2:3],e[,2:5])
        
        a <- be
        
        sds <- apply(a,2,sd)
        means <- apply(a,2,mean)
        
        er <- qnorm(0.975)*(sds/sqrt(nrow(a)))
        lw <- means - er
        hg <- means + er
        
        if (i==1) {
           plot(means[1:6],6:1,pch=19,col=colzF[i],yaxt='n',cex=1.8,cex.axis=1.2,ylab='',xlim=c(-1,1),cex.lab=1.2,xlab="Bias (Final Model - CV Model)")
        } else {
                points(means[1:6],6:1,pch=19,col=colzF[i],yaxt='n',cex=1.8,cex.axis=1.2,ylab='')
          }
        
        arrows(lw[1:6],6:1,hg[1:6],6:1,
               code=3,length=0,angle=90,lwd=2,col=colzF[i])
        abline(v=0,col='grey')
}

# miscalibration (total error) absolute value of bias

for (i in 1:3) {
  b <- diff.b.a[[i]]
  e <- diff.e.a[[i]]
  be <- cbind(b[,2:3],e[,2:5])
  
  a <- be
  
  sds <- apply(a,2,sd)
  means <- apply(a,2,mean)
  
  er <- qnorm(0.975)*(sds/sqrt(nrow(a)))
  lw <- means - er
  hg <- means + er
  
  if (i==1) {
    plot(means[1:6],6:1,pch=19,col=colzF[i],yaxt='n',cex=1.8,cex.axis=1.2,ylab='',xlim=c(0,1.3),cex.lab=1.2, xlab="Miscalibration")
    legend("bottomright",       
           legend = c("Random", "Block - Longitude","Block - Latitude"),
           col = c('azure3','darkkhaki','darkolivegreen'),
           pch = 19,cex=1.2)
    
  } else {
    points(means[1:6],6:1,pch=19,col=colzF[i],yaxt='n',cex=1.8,cex.axis=1.2,ylab='')
  }
  
  arrows(lw[1:6],6:1,hg[1:6],6:1,
         code=3,length=0,angle=90,lwd=2,col=colzF[i])
}

dev.off()

#----------------------------------



#----------------------------------------

#example of average species responses

#-------------------------------------------------
#change in probability with a change in prey

#average group-level parameters
a<- colMeans(meanparams.gold[grep('A',rownames(meanparams.gold)),])
b<- colMeans(meanparams.gold[grep('B',rownames(meanparams.gold)),])
m<- colMeans(meanparams.gold[grep('M',rownames(meanparams.gold)),])
r<- colMeans(meanparams.gold[grep('R',rownames(meanparams.gold)),])

# find average species from each group to use as example
av.a <- meanparams.gold[which(round(meanparams.gold[,'prey'],0)==round(a[2],0)),] #A18
av.b <- meanparams.gold[which(round(meanparams.gold[,'prey'],0)==round(b[2],0)),] #B207
av.m <- meanparams.gold[which(round(meanparams.gold[,'prey'],0)==round(m[2],0)),] #M11
av.r <- meanparams.gold[which(round(meanparams.gold[,'prey'],0)==round(r[2],0)),] #R106

species <- c("A18","B207","M11","R106")
sp <- species[1]
taxa <-  paste(taxo[which(taxo$ID%in%sp),5],taxo[which(taxo$ID%in%sp),6],sep=" ")
for (i in 2:4) {
  sp <- species[i]
  taxa <- c(taxa,paste(taxo[which(taxo$ID%in%sp),5],taxo[which(taxo$ID%in%sp),6],sep=" "))
}

taxa[3] <- "Vulpes lagopus"




load("InputData/Climate5kmfootprintv2.Rdata")

## first, generate partial responses for 'hypothetical' average species for each group
#colnames(clim)[10] <- 'foot'
#s <- spdat[which(spdat$ID%in%rownames(meanparams.gold)),]
#groups <- c('A','B','M','R')
#source("CodeHypotheticalAverageSpecies.R)


###  second, now for a selection of actual species that have average prey values..

# run predictions
source("Functions/FunctionsCrossVali.R")
runIndModPredsExamples(sp.list=species,clim=clim,bio1="MW2Prey",bio2="MW2CompR25",type="R25",metawebName="MW2")


#plot response curves and 95% credible intervals------

#opaque
colz2 <- c("ivory2","ivory3","antiquewhite4")
#transparent
c1 <- rgb(col2rgb("ivory2")[1],col2rgb("ivory2")[2],col2rgb("ivory2")[3],max = 255, alpha = 150)
c2 <- rgb(col2rgb("ivory3")[1],col2rgb("ivory3")[2],col2rgb("ivory3")[3],max = 255, alpha = 150)
c3 <- rgb(col2rgb("antiquewhite4")[1],col2rgb("antiquewhite4")[2],col2rgb("antiquewhite4")[3],max = 255, alpha = 150)
colzTr <- c(c1,c2,c3)

#group colors
group.col.full <- c("#e7298a","#7570b3","#d95f02","#1b9e77")



  
#----
pdf(file="ExampleResponseCurves.pdf",width=10,height=10)
par(mfrow=c(4,3))


for (i in 1:length(species)) {
  sp <- species[i] 
  group <- gsub("[0-9]*","",sp)
  dat <- meanparams.gold[grep(group,rownames(meanparams.gold)),]
  means <- apply(dat,2,mean)
  sds <- apply(dat,2,sd)
  
  er <- qnorm(0.975)*(sds/sqrt(nrow(dat)))
  lw <- means - er
  hg <- means + er
  spmean <- dat[grep(paste(sp,"$",sep=''),rownames(dat)),]
  
  #plot means for all species in each group AND the mean param estimate of the example species from each group
  if(i<4) {
    par(mar=c(1,5,0,0))
    
    plot(1:6,means[2:7],pch=21,xaxt='n',cex=1.8,cex.axis=1.2,cex.lab=1.2,ylab="Mean Parameter Estimate",ylim=c(-2.5,3.5),col=group.col.full[i])
      
  } else { 
      par(mar=c(2,5,0,0)) 
    plot(1:6,means[2:7],pch=21,cex=1.8,cex.axis=1.2,xaxt='n',
         ylim=c(-2,3.5),cex.lab=1.2,ylab='Mean Parameter Estimate',
         col=group.col.full[i])
    axis(1,at=1:6,labels=c('Prey','Competitor','Precip.','Min. temp','Prec. seas.','Footprint'),las=1,cex.axis=1.2)
    }
    
   arrows(1:6,lw[2:7],1:6,hg[2:7],
           code=3,length=0,angle=90,lwd=2,col=group.col.full[i])
    abline(h=0,col='grey')
    
    points(c(1.1,2.1,3.1,4.1,5.1,6.1),spmean[2:7],pch=19,col='antiquewhite4',xaxt='n',cex=1.8,cex.axis=1.2,
             cex.lab=1.2,ylab='',ylab='Mean Effect (Final Model)')
    
    #preds for reponse curves

    load(file=paste("ExampleModels/ExamplePredictions.MW2.R25.SpSelect",sp,"Rdata",sep="."))
    load(file=paste("ExampleModels/ExampleBetas.MW2.R25.SpSelect",sp,"Rdata",sep="."))
         
    env.new$pr <- exp(env.new$preds) / (1 + exp(env.new$preds))
    env.new$pr.lw <- exp(env.new$q2p5) / (1 + exp(env.new$q2p5))
    env.new$pr.up <- exp(env.new$q97p5) / (1 + exp(env.new$q97p5))
    env.new$b6div <- env.new$b6/100
    
    #b6 min temp
    if(i<4) {
      par(mar=c(1,5,0,0))
      plot(env.new[1:100,'pr']~env.new[1:100,'b6div'],
           xlim=c(min(env.new$b6div),max(env.new$b6div)),cex=1.2,cex.lab=1.2,cex.axis=1.2,xaxt='n',ylab="Pr. Occurrence",
           type='l',lwd=2,col='white',ylim=c(0,1))
      
    } else { 
      par(mar=c(2,5,0,0))
      plot(env.new[301:400,'pr']~env.new[301:400,'b6div'],
           xlim=c(min(env.new$b6div),max(env.new$b6div)),cex.lab=1.2,cex=1.4,cex.axis=1.2,ylab="Pr. Occurrence",
           type='l',lwd=2,col='white',ylim=c(0,1),xlab="Min.Temperature")
    }
   
    

    polygon(c(env.new[201:300,'b6div'],rev(env.new[201:300,'b6div'])), 
            c(env.new[201:300,'pr.up'],rev(env.new[201:300,'pr.lw'])),border=NA,col=colzTr[3])
    if (max(env.new$MW2Prey)>1) {
      polygon(c(env.new[101:200,'b6div'],rev(env.new[101:200,'b6div'])), 
              c(env.new[101:200,'pr.up'],rev(env.new[101:200,'pr.lw'])),border=NA,col=colzTr[2])
    }
    polygon(c(env.new[1:100,'b6div'],rev(env.new[1:100,'b6div'])), 
            c(env.new[1:100,'pr.up'],rev(env.new[1:100,'pr.lw'])),border=NA,col=colzTr[1])
    
    lines(env.new[1:100,'pr']~env.new[1:100,'b6div'],type='l',lwd=3,col=colz2[1])
    if (max(env.new$MW2Prey)>1) {
      lines(env.new[101:200,'pr']~env.new[101:200,'b6div'],type='l',lwd=3,col=colz2[2])
      }
    lines(env.new[201:300,'pr']~env.new[201:300,'b6div'],type='l',lwd=3,col=colz2[3])
    text(-2.1,0.1,paste('Max.Prey =',max(env.new$MW2Prey)),cex=1.4,col='grey35')
    
    #footprint
    if(i<4) {
      par(mar=c(1,2,0,3))
      plot(env.new[301:400,'pr']~env.new[301:400,'foot'],
           xlim=c(min(env.new$foot),max(env.new$foot)),cex=1.2,cex.axis=1.2,yaxt='n',xaxt='n',
           type='l',lwd=2,col='white',ylim=c(0,1),ylab="",xlab="Footprint")
      
    } else { 
        par(mar=c(2,2,0,3))
      plot(env.new[301:400,'pr']~env.new[301:400,'foot'],
           xlim=c(min(env.new$foot),max(env.new$foot)),cex=1.4,cex.axis=1.2,yaxt='n',
           type='l',lwd=2,col='white',ylim=c(0,1),ylab="",xlab="Footprint")
      }
    

    
    polygon(c(env.new[501:600,'foot'],rev(env.new[501:600,'foot'])), 
            c(env.new[501:600,'pr.up'],rev(env.new[501:600,'pr.lw'])),border=NA,col=colzTr[3])
    if (max(env.new$MW2Prey)>1) {
      polygon(c(env.new[401:500,'foot'],rev(env.new[401:500,'foot'])), 
            c(env.new[401:500,'pr.up'],rev(env.new[401:500,'pr.lw'])),border=NA,col=colzTr[2])
      }
    polygon(c(env.new[301:400,'foot'],rev(env.new[301:400,'foot'])), 
            c(env.new[301:400,'pr.up'],rev(env.new[301:400,'pr.lw'])),border=NA,col=colzTr[1])
    
    lines(env.new[301:400,'pr']~env.new[301:400,'foot'],type='l',lwd=3,col=colz2[1])
    if (max(env.new$MW2Prey)>1) {
      lines(env.new[401:500,'pr']~env.new[401:500,'foot'],type='l',lwd=3,col=colz2[2])
      }
    lines(env.new[501:600,'pr']~env.new[501:600,'foot'],type='l',lwd=3,col=colz2[3])
    if (i==4) {
      legend("bottomright",legend=c("Few prey", "Median prey", "Many prey"), 
             col=colz2,cex=1.4,lty=1,bty='n')
    }

    text(30,0.9,taxa[i],font=3,cex=1.4,col='grey35')
    rm(env.new)
}
dev.off()

#-------------------------------

## third, generate partial responses for 'hypothetical' average species for each group 
#for the model with prey richness and species richness
files <- list.files(pattern='IndivModels/SingleSp6predsLOGspeciesrichness.MW2.R25')

modstat <- lapply(files,function (x) get(load(x)))
names(modstat) <- gsub(".Rdata","",gsub("SingleSp6predsLOGspeciesrichness.MW2.R25.SpSelect.","",files))

params.sr <- sapply(modstat, "[",1)
meanparams.sr <- t(as.data.frame(sapply(params.sr, function(x) x[,1])))
rownames(meanparams.sr) <- gsub('.statistics','',rownames(meanparams.sr))
colnames(meanparams.sr) <- c('alpha','prey','sr','b12','b6','b15','foot')
#colnames(clim)[10] <- 'foot'
#s <- spdat[which(spdat$ID%in%rownames(meanparams.gold)),]
#groups <- c('A','B','M','R')
#source("CodeHypotheticalAverageSpecies.R)
## parameter estimates (MW 2)










## summarize mean beta values for all runs for appendix
#----------------------------------------------------------
## comparison 1 - metawebs 

#r1=MW1,R25 
#r2=gold-standard betas- parameter estimates drawn from full models (MW 2, R25)
#r1=MW4,RA


runs <- list('SingleSp6preds.MW1.R25','SingleSp6preds.MW2.R25','SingleSp6preds.MW4.RA')

dat <- sapply(runs, function (x) sumparams(list.files(pattern=x)))

#plot
colzF=c('#8dd3c7','#bebada','#80b1d3')
labz=c("MW1","MW2","MW-Globi")


pdf(file="AppendixBetaVarMetawebs.pdf",width=7,height=7)
plot.betas(runs=dat,labs=labz,colz=colzF,xlim=c(-3,3))
dev.off()
rm(runs,dat)




## comparison 2 - between different Biotic Variables (different competition (MW2 - R25,RA,R90,logged))----

runs <- list('SingleSp6preds.MW2.RA','SingleSp6preds.MW2.R25','SingleSp6preds.MW2.R90','SingleSp6predsLOG.MW2.R25')
colzF=c('#8dd3c7','#bebada','#80b1d3','royal blue1')
labz=c("All","25%","90%","log")

dat <- sapply(runs, function (x) sumparams(list.files(pattern=x)))


pdf(file="AppendixBetaVarBioticMetrics.pdf",width=7,height=7)
plot.betas(runs=dat,labs=labz,colz=colzF,xlim=c(-3,3))
dev.off()
rm(runs,dat)


## comparison 3 -between cross-validation runs-----

# full models (all predictors) - Random, Block lat, block long
meanparams.a <- vector(mode="list", length=3)
names(meanparams.a) <- c('Random','Blong','Blat')

pts <- c('SingleSp6predsCVtrain.MW2.R25.SpSelect.',
         'SingleSp6predsCVtrainSpatialEast.MW2.R25.SpSelect.',
         'SingleSp6predsCVtrainSpatialNorth.MW2.R25.SpSelect.')

nam <- rownames(b[[1]])


for (i in 1:3) {
  s <- list.files(pattern=pts[i])
  sub <- gsub("SingleSp6predsCVtrain.MW2.R25.SpSelect.","",s)
  sub <- gsub("SingleSp6predsCVtrainSpatialEast.MW2.R25.SpSelect.","",sub)
  sub <- gsub("SingleSp6predsCVtrainSpatialNorth.MW2.R25.SpSelect.","",sub)
  sub <- gsub(".Rdata","",sub)
  s <- s[sub%in%nam]
  comb <- lapply(s,function(x) get(load(x)))
  params <- sapply(comb, "[",1)
  e <- t(as.data.frame(sapply(params, function(x) x[,1])))
  sp <- gsub(pts[i],"",s)
  sp <- gsub(".Rdata","",sp)
  rownames(e) <- sp
  colnames(e) <- c('alpha','prey','comp','b12','b6','b15','foot')
  meanparams.a[[i]] <- e
  rm(s,e,sp)
}

a <- lapply(seq_along(meanparams.a), function(i) {
  dat <- meanparams.a[[i]]
  subset(dat, row.names(dat) %in% indx)
})

## bio only (b) and env only (e) 

be1 <- cbind(b[[1]],e[[1]][,2:5])
be2 <- cbind(b[[2]],e[[2]][,2:5])
be3 <- cbind(b[[3]],e[[3]][,2:5])

dat <- list(a[[1]],a[[2]],a[[3]],be1,be2,be3)


colzF=c('wheat2','skyblue2','seagreen2','wheat3','skyblue3','seagreen3')       
labz=c("Full-random","Full-block-long","Full-block-lat","Part-random","Part-block-long","Part-Block-lat")


pdf(file="AppendixBetaVarCrossValidation.pdf",width=7,height=7)
plot.betas(runs=dat,labs=labz,colz=colzF,xlim=c(-3,3))
dev.off()


## comparison 4. Different taxonomic groups (biotic versus abiotic variables from cross validation with random sites )
a.all <- a[[1]][grep("A",rownames(a[[1]])),]
b.all <- a[[1]][grep("B",rownames(a[[1]])),]
m.all <- a[[1]][grep("M",rownames(a[[1]])),]
r.all <- a[[1]][grep("R",rownames(a[[1]])),]

a.part <- be1[grep("A",rownames(be1)),]
b.part <- be1[grep("B",rownames(be1)),]
m.part <- be1[grep("M",rownames(be1)),]
r.part <- be1[grep("R",rownames(be1)),]



BirdCol <- rgb(col2rgb("#7570b3")[1],col2rgb("#7570b3")[2],col2rgb("#7570b3")[3],max = 255, alpha = 90)
MamCol <- rgb(col2rgb("#d95f02")[1],col2rgb("#d95f02")[2],col2rgb("#d95f02")[3],max = 255, alpha = 90)
AmphCol <- rgb(col2rgb("#e7298a")[1],col2rgb("#e7298a")[2],col2rgb("#e7298a")[3],max = 255, alpha = 90)
ReptCol <- rgb(col2rgb("#1b9e77")[1],col2rgb("#1b9e77")[2],col2rgb("#1b9e77")[3],max = 255, alpha = 90)



colzF=c("#e7298a","#7570b3","#d95f02","#1b9e77",AmphCol,BirdCol,MamCol,ReptCol)   
labz=c("Full-Amph","Full-Bird","Full-Mamm","Full-Rept","Part-Amph","Part-Bird","Part-Mamm","Part-Rept")

dat <- list(a.all,b.all,m.all,r.all,a.part,b.part,m.part,r.part)


pdf(file="AppendixBetaVarTaxonGroups.pdf",width=7,height=7)
plot.betas(runs=dat,labs=labz,colz=colzF,xlim=c(-4,4))
dev.off()



## comparison 5 . GBIF amphibians (MW 2, R25) with one and two occurrences per grid cell compared to original dist. data (logged)


datoriglog <- sumparams(list.files(pattern='SingleSp6predsLOG.MW2.R25'))
datorig <- sumparams(list.files(pattern='SingleSp6preds.MW2.R25'))


datgbif1 <- sumparams(list.files(pattern='SingleSp6predsLOGGBIF.MW2.R25.SpSelect.'))


datgbif2 <- sumparams(list.files(pattern='SingleSp6predsLOGGBIF.MW2.R25.SpSelect.'))

datorig <- datorig[which(rownames(datorig)%in%rownames(datgbif1)),]
datoriglog <- datoriglog[which(rownames(datoriglog)%in%rownames(datgbif1)),]

dat <- list(datorig,datoriglog,datgbif1,datgbif2)


colzF=c('#8dd3c7','#bebada','#80b1d3','royal blue1')
labz=c("Orig.","Orig-Log","GBIF-1","GBIF-2")

setwd("~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/Analyses3Feb20")
pdf(file="AppendixGBIFamphs.pdf",width=7,height=7)

plot.betas(runs=dat,labs=labz,colz=colzF,xlim=c(-4.5,4.5))

dev.off()



## comparison 6 . Models with species richness rather than comp (MW 2, R25) 

datsr <- sumparams(list.files(pattern='SingleSp6predsLOGspeciesrichness.MW2.R25'))

colzF=c('#8dd3c7')
labz=c("SpeciesRichness")

setwd("~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/Analyses3Feb20")
pdf(file="AppendixSpeciesRichness.pdf",width=7,height=7)

    a <- datsr
    sds <- apply(a,2,sd)
    means <- apply(a,2,mean)
    
    er <- qnorm(0.975)*(sds/sqrt(nrow(a)))
    lw <- means - er
    hg <- means + er
    par(mar=c(5,7,2,1))
    
    plot(means[1:7],7:1,pch=19,col=colzF,yaxt='n',cex=1.8,cex.axis=1.2,xlim=c(-2.5,3),ylab='',cex.lab=1.2,xlab="Mean Parameter Est.")

    arrows(lw[1:7],7:1,hg[1:7],7:1,
           code=3,length=0,angle=90,lwd=2,col=colzF)
    abline(v=0,col='grey')
    axis(2,at=7:1,labels=c('intercept','Prey','SR','Precip.','Min. temp','Prec. seas.','Footprint'),las=1,cex.axis=1.2)
    
dev.off()






