
## MODEL FUNCTIONS-------------------------------------------------------------

runIndModBRTCrossVal <- function(sp.list,clim,bio1,bio2,type,metawebName){
  
  
  for (x in 1:length(sp.list)) {
    # for (x in 1:3) {
    print(paste("Start mod",sp.list[x],"at",Sys.time()))
    
    sp <- sp.list[x]
    
    load(paste("IndivModels/",sp,".siteIntersAllWebs.Rdata",sep=""))
    
    if (is.null(siteInters)) {print(paste("mod",sp.list[x],"failed to load"))}
    if (is.null(siteInters)) next
    
    t <- do.call("cbind",siteInters)
    
    bio <- t[,c("Occur",bio1,bio2)] 
    if (max(bio[,1])==0|max(bio[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio[,1])==0|max(bio[,2])==0) next
    
    datt <- merge(bio,clim[,c(3,6:10)],by.x=0,by.y="PageName",all=F)
    datt <- datt[,-1]
    colnames(datt)[8] <- 'foot'
    
    # split data in half------
    samp <- sample(nrow(datt),nrow(datt)/2,replace=FALSE)
    dat.train <- datt[samp,]
    dat.test <- datt[-samp,]
    
    m <- gbm.step(data=dat.train,gbm.x = 2:8,gbm.y = 1,family = "bernoulli",tree.complexity = 3,
                     learning.rate = 0.001,bag.fraction = 0.7)
    
    b <- predict.gbm(m, dat.test[,2:8], n.trees=m$gbm.call$best.trees, type="response")
    pr.test <- cbind(dat.test[,1],b)
    colnames(pr.test) <- c("test.occur","probs")
    
    saveRDS(pr.test,file=paste("BRTs/SingleSp6predsCVtest",metawebName,type,"SpSelect",sp,"rds",sep="."))
    saveRDS(summary(m),file=paste("BRTs/SingleSp6predsCVtrainParams",metawebName,type,"SpSelect",sp,"rds",sep="."))
    
    rm(m,b,pr.test)
    
    #  2. now model with only enviro vars------------
    
    env.clim.train <- dat.train[,c(1,4:8)]
    env.clim.test <- dat.test[,c("b4","b12","b6","b15","foot")]
    
    m <- gbm.step(data=env.clim.train,gbm.x = 2:ncol(env.clim.train),gbm.y = 1,family = "bernoulli",tree.complexity = 4,
                     learning.rate = 0.002,bag.fraction = 0.6)
    
    b <- predict.gbm(m, env.clim.test, n.trees=m$gbm.call$best.trees, type="response")
    pr.test <- cbind(dat.test[,1],b)
    colnames(pr.test) <- c("test.occur","probs")
    
    saveRDS(pr.test,file=paste("BRTs/SingleSp6predsCVtestEnvOnly",metawebName,type,"SpSelect",sp,"rds",sep="."))
    
    rm(m,b,pr.test)
    
    # 3. now only biotic predictors----------
    
    biot.train <- dat.train[,c(1:3)]
    biot.test <- dat.test[,c("MW2Prey","MW2CompR25")]
    
    m <- gbm.step(data=biot.train,gbm.x = 2:3,gbm.y = 1,family = "bernoulli",tree.complexity = 4,
                     learning.rate = 0.002,bag.fraction = 0.6)
    
    b <- predict.gbm(m, biot.test, n.trees=m$gbm.call$best.trees, type="response")
    pr.test <- cbind(dat.test[,1],b)
    colnames(pr.test) <- c("test.occur","probs")
    
    saveRDS(pr.test,file=paste("BRTs/SingleSp6predsCVtestBioticOnly",metawebName,type,"SpSelect",sp,"rds",sep="."))
    
    rm(m,b,pr.test,dat.test,dat.train)
    
  }
}





### NOW run cross validation on spatial blocks

runIndModBRTCrossValSpatial <- function(sp.list,clim,bio1,bio2,type,metawebName){
  
  
  for (x in 1:length(sp.list)) {
    # for (x in 1:3) {
    print(paste("Start mod",sp.list[x],"at",Sys.time()))
    
    sp <- sp.list[x]
    
    load(paste("IndivModels/",sp,".siteIntersAllWebs.Rdata",sep=""))
    
    if (is.null(siteInters)) {print(paste("mod",sp.list[x],"failed to load"))}
    if (is.null(siteInters)) next
    
    t <- do.call("cbind",siteInters)
    
    datt <- merge(t,clim[,c(3:10)],by.x=0,by.y="PageName",all=F)
    colnames(datt)[25] <- 'foot'
    
    bio <- datt[,c(bio1,bio2)] 
    
    if (max(bio[,1])==0|max(bio[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio[,1])==0|max(bio[,2])==0) next
    ##
    
    
    ## NORTH-SOUTH split----------------------------------------------------------
    # split data in half- train on south, test on north
    medi <- median(datt[datt$Occur==1,'y'])
    dat <- datt[which(datt$y < medi),]
    dat.test <- datt[which(datt$y > medi),]
    
    #dat train
    bio <- dat[,c(bio1,bio2)] 
    cli <- dat[,c('b12','b6','b15','foot')]
    
    if (max(bio[,1])==0|max(bio[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio[,1])==0|max(bio[,2])==0) next
    
    occupancy <- dat[,"Occur"]
    
    env <- cbind(occupancy,bio,cli)
    
    #dat test
    bio.t <- dat.test[,c(bio1,bio2)] 
    cli.t <- dat.test[,c('b12','b6','b15','foot')]
    
    if (max(bio.t[,1])==0|max(bio.t[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio.t[,1])==0|max(bio.t[,2])==0) next
    
    occupancy.t <- dat.test[,c("Occur")]
    
    env.t <- cbind(bio.t,cli.t)
    
    m <- gbm.step(data=env,gbm.x = 2:ncol(env),gbm.y = 1,family = "bernoulli",tree.complexity = 4,
                     learning.rate = 0.001,bag.fraction = 0.6)
    
    b <- predict.gbm(m, env.t, n.trees=m$gbm.call$best.trees, type="response")
    pr.test <- cbind(occupancy.t,b)
    colnames(pr.test) <- c("test.occur","probs")
    
    saveRDS(pr.test,file=paste("BRTs/SingleSp6predsCVtestSpatialNorth",metawebName,type,"SpSelect",sp,"rds",sep="."))
    saveRDS(summary(m),file=paste("BRTs/SingleSp6predsCVtrainParamsSpatialNorth",metawebName,type,"SpSelect",sp,"rds",sep="."))
    
    rm(m,pr.test)
    
    #  2. now model with only enviro vars------------
    
    env.clim <- env[,c("occupancy","b12","b6","b15","foot")]
    env.clim.t <- env.t[,c("b12","b6","b15","foot")]
    
    m <- gbm.step(data=env.clim,gbm.x = 2:ncol(env.clim),gbm.y = 1,family = "bernoulli",tree.complexity = 4,
                  learning.rate = 0.001,bag.fraction = 0.6)
    
    b <- predict.gbm(m, env.clim.t, n.trees=m$gbm.call$best.trees, type="response")
    pr.test <- cbind(occupancy.t,b)
    colnames(pr.test) <- c("test.occur","probs")
    
    saveRDS(pr.test,file=paste("BRTs/SingleSp6predsCVtestEnvOnlySpatialNorth",metawebName,type,"SpSelect",sp,"rds",sep="."))

    
    # 3. now only biotic predictors----------
    
    env.biot <- env[,c("occupancy","MW2Prey","MW2CompR25")]
    env.biot.t <- env.t[,c("MW2Prey","MW2CompR25")]
    
    m <- gbm.step(data=env.biot,gbm.x = 2:ncol(env.biot),gbm.y = 1,family = "bernoulli",tree.complexity = 4,
                  learning.rate = 0.001,bag.fraction = 0.6)
    
    b <- predict.gbm(m, env.biot.t, n.trees=m$gbm.call$best.trees, type="response")
    pr.test <- cbind(occupancy.t,b)
    colnames(pr.test) <- c("test.occur","probs")
    
    saveRDS(pr.test,file=paste("BRTs/SingleSp6predsCVtestBioticOnlySpatialNorth",metawebName,type,"SpSelect",sp,"rds",sep="."))
    
    rm(b,pr.test,m)
    
    
    ## EAST-WEST split----------------------------------------------------------
    
    
    # split data in half- train on west, test on east
    medi <- median(datt[datt$Occur==1,'x'])
    dat <- datt[which(datt$x < medi),]
    dat.test <- datt[which(datt$x > medi),]
    
    #dat train
    bio <- dat[,c(bio1,bio2)] 
    cli <- dat[,c('b12','b6','b15','foot')]
    
    if (max(bio[,1])==0|max(bio[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio[,1])==0|max(bio[,2])==0) next
    
    ##
    occupancy <- dat[,"Occur"]
    
    env <- cbind(occupancy,bio,cli)
    
    #dat test
    bio.t <- dat.test[,c(bio1,bio2)] 
    cli.t <- dat.test[,c('b12','b6','b15','foot')]
    
    if (max(bio.t[,1])==0|max(bio.t[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio.t[,1])==0|max(bio.t[,2])==0) next
    
    occupancy.t <- dat.test[,c("Occur")]
    
    env.t <- cbind(bio.t,cli.t)
    
    m <- gbm.step(data=env,gbm.x = 2:ncol(env),gbm.y = 1,family = "bernoulli",tree.complexity = 4,
                  learning.rate = 0.001,bag.fraction = 0.6)
    
    b <- predict.gbm(m, env.t, n.trees=m$gbm.call$best.trees, type="response")
    pr.test <- cbind(occupancy.t,b)
    colnames(pr.test) <- c("test.occur","probs")
    
    saveRDS(pr.test,file=paste("BRTs/SingleSp6predsCVtestSpatialEast",metawebName,type,"SpSelect",sp,"rds",sep="."))
    saveRDS(summary(m),file=paste("BRTs/SingleSp6predsCVtrainParamsSpatialEast",metawebName,type,"SpSelect",sp,"rds",sep="."))
    
    rm(m,pr.test)
    
    #  2. now model with only enviro vars------------
    
    env.clim <- env[,c("occupancy","b12","b6","b15","foot")]
    env.clim.t <- env.t[,c("b12","b6","b15","foot")]
    
    m <- gbm.step(data=env.clim,gbm.x = 2:ncol(env.clim),gbm.y = 1,family = "bernoulli",tree.complexity = 4,
                  learning.rate = 0.001,bag.fraction = 0.6)
    
    b <- predict.gbm(m, env.clim.t, n.trees=m$gbm.call$best.trees, type="response")
    pr.test <- cbind(occupancy.t,b)
    colnames(pr.test) <- c("test.occur","probs")
    
    saveRDS(pr.test,file=paste("BRTs/SingleSp6predsCVtestEnvOnlySpatialEast",metawebName,type,"SpSelect",sp,"rds",sep="."))
    
    
    # 3. now only biotic predictors----------
    
    env.biot <- env[,c("occupancy","MW2Prey","MW2CompR25")]
    env.biot.t <- env.t[,c("MW2Prey","MW2CompR25")]
    
    m <- gbm.step(data=env.biot,gbm.x = 2:ncol(env.biot),gbm.y = 1,family = "bernoulli",tree.complexity = 4,
                  learning.rate = 0.001,bag.fraction = 0.6)
    
    b <- predict.gbm(m, env.biot.t, n.trees=m$gbm.call$best.trees, type="response")
    pr.test <- cbind(occupancy.t,b)
    colnames(pr.test) <- c("test.occur","probs")
    
    saveRDS(pr.test,file=paste("BRTs/SingleSp6predsCVtestBioticOnlySpatialEast",metawebName,type,"SpSelect",sp,"rds",sep="."))
    
    rm(b,pr.test,m)
    
    
    rm(m,b,pr.test,dat,dat.test,occupancy,occupancy.t)
  }
}
