
library(Matrix)
library(reshape2)


# prep

#1  SET 1 - ANY SPECIES THAT HAS PREY 
setwd("~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/Analyses3Feb20/InputData/")

##bio4TempSeas, bio12AnPrecip
load("Climate5kmfootprintv2.Rdata")
colnames(clim)[10] <- "foot"

load("SiteBySpecies21Feb.Rdata")
rownames(clim) <- clim$PageName
site.by.species <- site.by.species[which(rownames(site.by.species) %in% rownames(clim)),]

# first for any species with prey in the most potential metaweb
load("pred.by.prey.Metaweb1.Rdata")
predator.list <- apply(pred.by.prey,1,sum).
predator.list <- names(predator.list[which(predator.list>0)]) #367 non basal species


# 5 runs for all species that have prey
#MW1 - R25, ACR25
#MW2 - R25 (default), R90,
#MW4globi - RA
MW1 <- pred.by.prey
rm(pred.by.prey)
load("pred.by.prey.Metaweb2Realized.Rdata")
MW2 <- pred.by.prey
rm(pred.by.prey)
load("pred.by.prey.Metaweb4Globi.Rdata")
MW4 <- pred.by.prey
rm(pred.by.prey)


p <- predator.list


source("FunctionsPythonModels.R")

## run for all species
runBioticInputsIndSpeciesWPrey(MW1=MW1,MW2=MW2,MW4=MW4,site.by.species=site.by.species,sp.list=p)

## run for species to predict to entire map (warning: very large operation)
runBioticInputsMaps(sp.list=sp.list,site.by.species=site.by.species,clim=clim,metaweb=pred.by.prey)


# run models------------------------------

reticulate::use_condaenv("r-tensorflow", required = TRUE)


##### Load packages and data ####


library(greta)
library(igraph)
library(DiagrammeR)
library(bayesplot)
library(Matrix)
library(reshape2)

setwd("~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/Analyses21Feb19")
source("FunctionsPythonModels.R")


#climate inputs

load("Climate5kmfootprintv2.Rdata")



#---------
load("pred.by.prey.Metaweb4Globi.Rdata")
sp.list <- names(rowSums(pred.by.prey)[which(rowSums(pred.by.prey)>0)])



# dir for individual mods
setwd("~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/Analyses3Feb20/IndivModels/")


## merge with climate data
sp.list <- gsub('.siteIntersAllWebs.Rdata','',list.files(pattern="siteIntersAll")) #367 


d <- gsub('SingleSp6preds.MW4.RA.SpSelect.','',list.files(pattern="SingleSp6preds.MW4.")) #367 
d <- gsub('.Rdata','',d) #367 

d <- gsub('SingleSp6preds.MW1.R25.SpSelect.','',list.files(pattern="SingleSp6preds.MW1.")) #367 
d <- gsub('.Rdata','',d) #367 

d <- gsub('SingleSp6preds.MW4.ACompRA.SpSelect.','',list.files(pattern="SingleSp6preds.MW4.ACompRA")) #367 
d <- gsub('.Rdata','',d) #367 

sp.list <- sp.list[!(sp.list %in% d)]
sp.list <- sp.list[!(sp.list=="B227")]

sp.list <- sp.list[-1]


sp.list <- sp.list[42:length(sp.list)]





# dir for individual mods
setwd("~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/Analyses3Feb20/IndivModels")

# RUNS MODELS
runIndModPython(sp.list=sp.list,clim=clim,bio1="MW4Prey",bio2="MW4ACompRA",type="ACompRA",metawebName="MW4")

runIndModPython(sp.list=sp.list,clim=clim,bio1="MW1Prey",bio2="MW1CompR25",type="R25",metawebName="MW1")

runIndModPython(sp.list=sp.list,clim=clim,bio1="MW2Prey",bio2="MW2CompRA",type="RA",metawebName="MW2")

runIndModPython(sp.list=sp.list,clim=clim,bio1="MW2Prey",bio2="MW2CompR25",type="R25",metawebName="MW2")

runIndModPython(sp.list=sp.list,clim=clim,bio1="MW4Prey",bio2="MW4CompRA",type="RA",metawebName="MW4")



# version with prey and competitor number logged
runIndModPythonLOGpreycomp(sp.list=sp.list,clim=clim,bio1="MW2Prey",bio2="MW2CompR25",type="R25",metawebName="MW2")

runIndModPythonLOGpreycomp(sp.list=sp.list,clim=clim,bio1="MW2Prey",bio2="MW2CompR25",type="R25",metawebName="MW2")


# version with prey and species richness rather than competitor richness
setwd("~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/Analyses3Feb20/InputData")
load("SiteBySpecies21Feb.Rdata")
sr <- rowSums(site.by.species)


runIndModPythonLOGpreySR(sp.list=sp.list,clim=clim,bio1="MW2Prey",richness=sr,type="R25",metawebName="MW2")


#---------------------------------------------------------------------
## make predictions for map
# dir for individual mods
setwd("~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/Analyses3Feb20/PredictMaps/")
sp.list <- gsub('.siteIntersMap.Rdata','',list.files(pattern="siteIntersMap")) #367 


setwd("~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/Analyses3Feb20")
runIndModPredsMaps(sp.list,clim)
#----------------------------------------------

#create species richness predictions for 'as-is', plus 25% prey and plut 50% prey in cells

filz <- list.files(pattern="Predictions.")
#predsActual <- sapply(filz, function(x) get(load(x))[,'pr.map']  )
#memory exhausted..

prmap <- get(load(filz[1]))[,c('pr.map','pr.map25','pr.map50')]

for (i in 2:length(filz)) {
   f <- get(load(filz[i]))[,c('pr.map','pr.map25','pr.map50')]
   prmap <- prmap + f
   rm(f)
}

nm <- get(load(filz[1]))[,c('Row.names')]
rownames(prmap) <- nm
prmap <- apply(prmap,2,function(x) round(x,digits=4))


# input layers for spatial grid and lakes
load('~/Dropbox/Manuscripts/EuropeanFoodWebCons/manuscript/Manuscript/ScienceNewAnalyses/MapPrepNewApr18.Rdata')

grid <- raster("InputData/reference_grid_5km.img")

vals <- merge(all,prmap,by.x="PageName",by.y=0,all.x=T)
valsord <- vals[order(vals[,'uniqueID.x']),]


#SR with lake mask
valsord$srmask <- ifelse(!is.na(valsord$lakes),NA,valsord$pr.map)
valsord$srmask <- ifelse(is.na(valsord$studyextent),NA,valsord$srmask)
#increase prey rich with mask
valsord$sr25mask <- ifelse(!is.na(valsord$lakes),NA,valsord$pr.map25)
valsord$sr25mask <- ifelse(is.na(valsord$studyextent),NA,valsord$sr25mask)
#increase prey rich with mask by 50%
valsord$sr50mask <- ifelse(!is.na(valsord$lakes),NA,valsord$pr.map50)
valsord$sr50mask <- ifelse(is.na(valsord$studyextent),NA,valsord$sr50mask)




# percent increase in richness

prActual <- grid
prActual[] <- valsord[,'srmask']

pr25 <- grid
pr25[] <- valsord[,"sr25mask"]

diff <- pr25 - prActual
diffperc <- (diff/prActual)*100



pdf(file="FigureMapsSR.pdf")
par(mfrow=c(1,1))

breakpoints <- seq(0,250,length.out=101)
colors <- colorRampPalette(brewer.pal(9,'PuBuGn'))(100)
plot(prActual,breaks=breakpoints,col=colors,axes=F,box=F,legend=T,main="SR")
#plot(pr25,breaks=breakpoints,col=colors,axes=F,box=F,legend=T,main="SR - 25% more prey")
dev.off()

pdf(file="FigureMapsSRplus25prey.pdf")
colors <- colorRampPalette(brewer.pal(9,'YlGnBu'))(100)
breakpoints <- seq(0,60,length.out=101)
plot(diffperc,breaks=breakpoints,col=colors,axes=F,box=F,legend=T,main="SR - 25% more prey")
dev.off()

vals$diff <- vals$pr.map25 - vals$pr.map
vals$diffperc <- (vals$diff/vals$pr.map)*100
samp <- vals[!is.na(vals$pr.map25),]
samp <- samp[sample(nrow(samp),100000),]



#lat <- clim[,c("PageName","x")]
#valsm <- merge(vals,lat,by="PageName")

pdf(file="FigureSRPreyIncreaseScatter.pdf")

colors <- colorRampPalette(brewer.pal(9,'PuBuGn'))(100)
breakpoints <- seq(0,60,length.out=101)
# Rank variable for colour assignment
samp$order = findInterval(samp$diffperc, breakpoints)

col2rgb("#519DC8")

trans <- rgb(81, 157, 200, max = 255, alpha = 7)
par(mar=c(4,6,4,4))
plot(samp$diffperc~samp$pr.map,col=trans,pch=19,cex.axis=1.8,cex.lab=1.8,ylab="% Increase in SR with 25% prey increase",xlab="SR")

dev.off()
##

pdf(file="FigureSRScatterAppendix.pdf")
boxplot(valsord$srmask,valsord$sr25mask,valsord$sr50mask,na.rm=T,ylab='Species Richness')
axis(side=1,at=c(1,2,3),labels=c('Current','Increase 25%','Increase 50%'))
dev.off()





#---
datcors <- vector("list", length(sp.list))
names(datcors) <- sp.list

for (x in 1:length(sp.list)) {
  
  sp <- sp.list[x]

  load(paste(sp,"siteIntersAllWebs.Rdata",sep="."))
  
  
   t <- do.call("cbind",siteInters)
   
   dat <- merge(t,clim[,c(3,6:10)],by.x=0,by.y="PageName",all=F)
   colnames(dat)[23] <- 'foot'
   
   datcors[[x]] <- cor(dat[,c(7:12,19:ncol(dat))])

   }
   
 
add <- function(x) Reduce("mean", x)
add(datcors)

sumcor <- Reduce('mean',datcors)

cadd <- function(x) Reduce("+", x, accumulate = TRUE)


#--------''

   



#----

runbioticPythoncomp(bioclim,comptype="compRA",metawebName='MetaWeb4SiteSelect')
runbioticPythoncomp(bioclim,comptype="compIn",metawebName='MetaWeb4SiteSelect')
runbioticPythoncomp(bioclim,comptype="compACra",metawebName='MetaWeb4SpSelect')







.



