

## functions used for cross-validation

sumparams <- function(files) {
  modstat <- lapply(files,function (x) get(load(x)))
  names(modstat) <- gsub(".Rdata","",gsub("SingleSp6preds.MW2.R25.SpSelect.","",files))
  
  params <- sapply(modstat, "[",1)
  meanparams <- t(as.data.frame(sapply(params, function(x) x[,1])))
  rownames(meanparams) <- gsub('.statistics','',rownames(meanparams))
  rownames(meanparams) <- gsub('SingleSp6predsLOGGBIF.MW2.R25.SpSelect.','',rownames(meanparams))
  rownames(meanparams) <- gsub('SingleSp6predsLOG.MW2.R25.SpSelect.','',rownames(meanparams))
  rownames(meanparams) <- gsub('SingleSp6predsLOGspeciesrichness.MW2.R25.SpSelect.','',rownames(meanparams))
  colnames(meanparams) <- c('alpha','prey','comp','b12','b6','b15','foot')
  return(meanparams)
}



plot.betas <- function(runs,labs,colz,xlim) {
  
  for (i in 1:length(runs)) {
    a <- runs[[i]]
    sds <- apply(a,2,sd)
    means <- apply(a,2,mean)
    
    er <- qnorm(0.975)*(sds/sqrt(nrow(a)))
    lw <- means - er
    hg <- means + er
    par(mar=c(5,7,2,1))
    
    if (i==1){
      plot(means[1:7],7:1,pch=19,col=colzF[i],yaxt='n',cex=1.8,cex.axis=1.2,ylab='',xlim=xlim,cex.lab=1.2,xlab="Mean Parameter Est.")
    } else {
      points(means[1:7],7:1,pch=19,col=colzF[i],yaxt='n',cex=1.8,cex.axis=1.2,ylab='')
    }
    
    arrows(lw[1:7],7:1,hg[1:7],7:1,
           code=3,length=0,angle=90,lwd=2,col=colzF[i])
    abline(v=0,col='grey')
    axis(2,at=7:1,labels=c('intercept','Prey','Competitor','Precip.','Min. temp','Prec. seas.','Footprint'),las=1,cex.axis=1.2)
    
  }
  par(mar=c(1,1,1,1))
  legend(1.3,2.5,       
         legend = labs,
         col = colzF,
         pch = 19,cex=1.1)
}


pr.fn <- function (i) exp(i) / (1 + exp(i))


combin <- function(x) {
  data.frame(
    meanlogit = x[['statistics']][,'Mean'],
    occur = x[['occur']]
  )}

agg <- function(prefx,n) {
    filz <- list.files(pattern=prefx)
    sp <- gsub(paste(prefx,'.R25.SpSelect.',sep=""),'',filz)
    sp <- gsub('.Rdata','',sp)
    
    comb <- lapply(filz,function(x) get(load(x)))
    names(comb) <- sp
    l <- lapply(comb, FUN=combin)
    logit <- do.call("rbind",l)
    logit$pr <- pr.fn(logit$meanlogit)
    return(logit)
  }



#make prediction object to use functions from rocr - must make separate list for any new changes in occur
## and for each species separately..

rocob <- function(splong,occur,pr,metric){
  auc <- data.frame(matrix(NA,nrow=length(unique(splong)),ncol=9))
  rownames(auc) <- unique(splong)
  
  for (i in 1:length(unique(splong))) {
    species <- unique(splong)[i]
    rows <- which(splong==species)
    O <- sapply(occur,'[',rows)
    if (min(colSums(O))==0) next
    P <- sapply(pr,'[',rows)
    
    pred.ob <- prediction(P,O)
    perf <- performance(pred.ob,metric)
    auc[species,] <- rbind(perf@y.values)
    rm(species,O,P,pred.ob,perf)
  }
  return(auc)
}

rocobBRT <- function(prefx,path=pth,metric,sp){
  
  auc <- data.frame(matrix(NA,nrow=length(sp),ncol=9))
  rownames(auc) <- sp
  
  # different CV runs
  for (i in 1:length(prefx)) {
    filz <- list.files(path=pth,pattern=prefx[i])  
    species <- gsub(paste(prefx[i],'.R25.SpSelect.',sep=""),'',filz)
    species <- gsub('.Rdata','',species)
    
    for (x in 1:length(species)){
      if (i==2|i==3) {
        load(paste(pth,filz[x],sep="/"))
        dat <- pr.test
      } else { 
        dat <- readRDS(paste(pth,filz[x],sep="/")) 
      }
      
      if (length(unique(dat[,1]))!=2) next
      
      pred.ob <- prediction(dat[,2],dat[,1])
      perf <- performance(pred.ob,metric)
      auc[species==rownames(auc)[x],i] <- perf@y.values
    }  
  }
  return(auc)   
}




#-----
plot.preds.auc <- function(nmz,x,m,colz,colzF,x2,ylab) {
  
  plot(x,m[,nmz[1]],col='grey',pch=19,ylab=ylab,xlab="log(Range Size)",ylim=c(0,1))
  
  #plot lines
  for(i in 1:length(nmz)){
    y=nmz[i]
    points(x,m[,y],col=colz[i],pch=19)
  }
  for(i in 1:length(nmz)){
    y=nmz[i]
    mod <- glm(m[,y]~x + x2 ,family=stats::binomial)
    if (summary(mod)$coef[2,4]<0.05) {
      p <- predict(mod,pred.mat, type='response')
      lines(seq(min(x),max(x),length.out=100),p,col=colzF[i],lwd=3)
    } 
    if (summary(mod)$coef[2,4]>=0.05) {
      p <- predict(mod,pred.mat, type='response')
      lines(seq(min(x),max(x),length.out=100),p,col=colzF[i],lwd=3,lty=2)
    } 
  }
  
}


## MODEL FUNCTIONS-------------------------------------------------------------

runIndModPythonCrossVal <- function(sp.list,clim,bio1,bio2,type,metawebName){
  
  
  for (x in 1:length(sp.list)) {
    # for (x in 1:3) {
    print(paste("Start mod",sp.list[x],"at",Sys.time()))
    
    sp <- sp.list[x]
    
    load(paste("IndivModels/",sp,".siteIntersAllWebs.Rdata",sep=""))
    
    if (is.null(siteInters)) {print(paste("mod",sp.list[x],"failed to load"))}
    if (is.null(siteInters)) next
    
    t <- do.call("cbind",siteInters)
    
    datt <- merge(t,clim[,c(3,6:10)],by.x=0,by.y="PageName",all=F)
    colnames(datt)[23] <- 'foot'
    
    bio <- datt[,c(bio1,bio2)] 
    
    if (max(bio[,1])==0|max(bio[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio[,1])==0|max(bio[,2])==0) next
    ##
    
    # split data in half
    samp <- sample(nrow(datt),nrow(datt)/2,replace=FALSE)
    dat <- datt[samp,]
    dat.test <- datt[-samp,]
    
    #dat train
    bio <- dat[,c(bio1,bio2)] 
    cli <- dat[,c('b12','b6','b15','foot')]
    
    ##
    occupancy <- dat[,"Occur"]
    
    env <- cbind(bio,cli)
    env <- apply(env,2,scale)
    
    #dat test
    bio.t <- dat.test[,c(bio1,bio2)] 
    cli.t <- dat.test[,c('b12','b6','b15','foot')]
    
    occupancy.t <- dat.test[,c("Occur")]
    
    env.t <- cbind(bio.t,cli.t)
    env.t <- apply(env.t,2,scale)
    
    ## model
    n_env <- ncol(env)
    n_sites <- nrow(env)
    
    #greta
    alpha <- normal(0, 1)
    beta <- normal(0, 1, dim = n_env)
    
    # logit-linear model
    linear_predictor <- alpha + env %*% beta
    p <- ilogit(linear_predictor)
    
    # distribution (likelihood) over observed values
    distribution(occupancy) <- bernoulli(p)
    
    #defining model
    m <- model(alpha,beta)
    
    #sampling 10.22
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 5000)
    
    sumdraws <- summary(draws)
    
    save(sumdraws,file=paste("CrossValidationModels/SingleSp6predsCVtrain",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    New <- alpha + env.t %*% beta
    New.draws <- calculate(New,draws)
    
    logi <- summary(New.draws)
    logi$occur <- occupancy.t
    
    save(logi,file=paste("CrossValidationModels/SingleSp6predsCVtest",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(alpha,beta,draws,sumdraws,New,New.draws,logi)
    
    #  2. now model with only enviro vars------------
    
    env.clim <- env[,c("b12","b6","b15","foot")]
    env.clim.t <- env.t[,c("b12","b6","b15","foot")]
    
    ## model
    n_env <- ncol(env.clim)
    n_sites <- nrow(env.clim)
    
    #greta
    alpha <- normal(0, 1)
    beta <- normal(0, 1, dim = n_env)
    
    # logit-linear model
    linear_predictor <- alpha + env.clim %*% beta
    p <- ilogit(linear_predictor)
    
    # distribution (likelihood) over observed values
    distribution(occupancy) <- bernoulli(p)
    
    #defining model
    m <- model(alpha,beta)
    
    #sampling 10.22
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 5000)
    
    sumdraws <- summary(draws)
    
    save(sumdraws,file=paste("CrossValidationModels/SingleSp6predsCVtrainEnvOnly",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    New <- alpha + env.clim.t %*% beta
    New.draws <- calculate(New,draws)
    
    logi <- summary(New.draws)
    logi$occur <- occupancy.t
    
    save(logi,file=paste("CrossValidationModels/SingleSp6predsCVtestEnvOnly",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(alpha,beta,draws,sumdraws,New,New.draws,logi,env.clim,env.clim.t)
    
    # 3. now only biotic predictors----------
    
    env.biot <- env[,c("MW2Prey","MW2CompR25")]
    env.biot.t <- env.t[,c("MW2Prey","MW2CompR25")]
    
    ## model
    n_env <- ncol(env.biot)
    n_sites <- nrow(env.biot)
    
    #greta
    alpha <- normal(0, 1)
    beta <- normal(0, 1, dim = n_env)
    
    # logit-linear model
    linear_predictor <- alpha + env.biot %*% beta
    p <- ilogit(linear_predictor)
    
    # distribution (likelihood) over observed values
    distribution(occupancy) <- bernoulli(p)
    
    #defining model
    m <- model(alpha,beta)
    
    #sampling 10.22
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 5000)
    
    sumdraws <- summary(draws)
    
    save(sumdraws,file=paste("CrossValidationModels/SingleSp6predsCVtrainBioticOnly",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    New <- alpha + env.biot.t %*% beta
    New.draws <- calculate(New,draws)
    
    logi <- summary(New.draws)
    logi$occur <- occupancy.t
    
    save(logi,file=paste("CrossValidationModels/SingleSp6predsCVtestBioticOnly",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(alpha,beta,draws,sumdraws,New,New.draws,logi,env,env.biot)
  }
}





### NOW run cross validation on spatial blocks

runIndModPythonCrossValSpatial <- function(sp.list,clim,bio1,bio2,type,metawebName){
  
  
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
    
    ##
    occupancy <- dat[,"Occur"]
    
    env <- cbind(bio,cli)
    env <- apply(env,2,scale)
    
    #dat test
    bio.t <- dat.test[,c(bio1,bio2)] 
    cli.t <- dat.test[,c('b12','b6','b15','foot')]
    
    if (max(bio.t[,1])==0|max(bio.t[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio.t[,1])==0|max(bio.t[,2])==0) next
    
    occupancy.t <- dat.test[,c("Occur")]
    
    env.t <- cbind(bio.t,cli.t)
    env.t <- apply(env.t,2,scale)
    
    ## model
    n_env <- ncol(env)
    n_sites <- nrow(env)
    
    #greta
    alpha <- normal(0, 1)
    beta <- normal(0, 1, dim = n_env)
    
    # logit-linear model
    linear_predictor <- alpha + env %*% beta
    p <- ilogit(linear_predictor)
    
    # distribution (likelihood) over observed values
    distribution(occupancy) <- bernoulli(p)
    
    #defining model
    m <- model(alpha,beta)
    
    #sampling 10.22
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 5000)
    
    sumdraws <- summary(draws)
    
    save(sumdraws,file=paste("CrossValidationModels/SingleSp6predsCVtrainSpatialNorth",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    New <- alpha + env.t %*% beta
    New.draws <- calculate(New,draws)
    
    logi <- summary(New.draws)
    logi$occur <- occupancy.t
    
    save(logi,file=paste("CrossValidationModels/SingleSp6predsCVtestSpatialNorth",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(alpha,beta,draws,sumdraws,New,New.draws,logi)
    
    #  2. now model with only enviro vars------------
    
    env.clim <- env[,c("b12","b6","b15","foot")]
    env.clim.t <- env.t[,c("b12","b6","b15","foot")]
    
    ## model
    n_env <- ncol(env.clim)
    n_sites <- nrow(env.clim)
    
    #greta
    alpha <- normal(0, 1)
    beta <- normal(0, 1, dim = n_env)
    
    # logit-linear model
    linear_predictor <- alpha + env.clim %*% beta
    p <- ilogit(linear_predictor)
    
    # distribution (likelihood) over observed values
    distribution(occupancy) <- bernoulli(p)
    
    #defining model
    m <- model(alpha,beta)
    
    #sampling 10.22
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 5000)
    
    sumdraws <- summary(draws)
    
    save(sumdraws,file=paste("CrossValidationModels/SingleSp6predsCVtrainEnvOnlySpatialNorth",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    New <- alpha + env.clim.t %*% beta
    New.draws <- calculate(New,draws)
    
    logi <- summary(New.draws)
    logi$occur <- occupancy.t
    
    save(logi,file=paste("CrossValidationModels/SingleSp6predsCVtestEnvOnlySpatialNorth",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(alpha,beta,draws,sumdraws,New,New.draws,logi,env.clim,env.clim.t)
    
    # 3. now only biotic predictors----------
    
    env.biot <- env[,c("MW2Prey","MW2CompR25")]
    env.biot.t <- env.t[,c("MW2Prey","MW2CompR25")]
    
    ## model
    n_env <- ncol(env.biot)
    n_sites <- nrow(env.biot)
    
    #greta
    alpha <- normal(0, 1)
    beta <- normal(0, 1, dim = n_env)
    
    # logit-linear model
    linear_predictor <- alpha + env.biot %*% beta
    p <- ilogit(linear_predictor)
    
    # distribution (likelihood) over observed values
    distribution(occupancy) <- bernoulli(p)
    
    #defining model
    m <- model(alpha,beta)
    
    #sampling 10.22
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 5000)
    
    sumdraws <- summary(draws)
    
    save(sumdraws,file=paste("CrossValidationModels/SingleSp6predsCVtrainBioticOnlySpatialNorth",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    New <- alpha + env.biot.t %*% beta
    New.draws <- calculate(New,draws)
    
    logi <- summary(New.draws)
    logi$occur <- occupancy.t
    
    save(logi,file=paste("CrossValidationModels/SingleSp6predsCVtestBioticOnlySpatialNorth",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(medi,dat,dat.test)
    #---------
    
    
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
    
    env <- cbind(bio,cli)
    env <- apply(env,2,scale)
    
    #dat test
    bio.t <- dat.test[,c(bio1,bio2)] 
    cli.t <- dat.test[,c('b12','b6','b15','foot')]
    
    if (max(bio.t[,1])==0|max(bio.t[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio.t[,1])==0|max(bio.t[,2])==0) next
    
    occupancy.t <- dat.test[,c("Occur")]
    
    env.t <- cbind(bio.t,cli.t)
    env.t <- apply(env.t,2,scale)
    
    ## model
    n_env <- ncol(env)
    n_sites <- nrow(env)
    
    #greta
    alpha <- normal(0, 1)
    beta <- normal(0, 1, dim = n_env)
    
    # logit-linear model
    linear_predictor <- alpha + env %*% beta
    p <- ilogit(linear_predictor)
    
    # distribution (likelihood) over observed values
    distribution(occupancy) <- bernoulli(p)
    
    #defining model
    m <- model(alpha,beta)
    
    #sampling 10.22
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 5000)
    
    sumdraws <- summary(draws)
    
    save(sumdraws,file=paste("CrossValidationModels/SingleSp6predsCVtrainSpatialEast",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    New <- alpha + env.t %*% beta
    New.draws <- calculate(New,draws)
    
    logi <- summary(New.draws)
    logi$occur <- occupancy.t
    
    save(logi,file=paste("CrossValidationModels/SingleSp6predsCVtestSpatialEast",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(alpha,beta,draws,sumdraws,New,New.draws,logi)
    
    #  2. now model with only enviro vars------------
    
    env.clim <- env[,c("b12","b6","b15","foot")]
    env.clim.t <- env.t[,c("b12","b6","b15","foot")]
    
    ## model
    n_env <- ncol(env.clim)
    n_sites <- nrow(env.clim)
    
    #greta
    alpha <- normal(0, 1)
    beta <- normal(0, 1, dim = n_env)
    
    # logit-linear model
    linear_predictor <- alpha + env.clim %*% beta
    p <- ilogit(linear_predictor)
    
    # distribution (likelihood) over observed values
    distribution(occupancy) <- bernoulli(p)
    
    #defining model
    m <- model(alpha,beta)
    
    #sampling 10.22
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 5000)
    
    sumdraws <- summary(draws)
    
    save(sumdraws,file=paste("CrossValidationModels/SingleSp6predsCVtrainEnvOnlySpatialEast",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    New <- alpha + env.clim.t %*% beta
    New.draws <- calculate(New,draws)
    
    logi <- summary(New.draws)
    logi$occur <- occupancy.t
    
    save(logi,file=paste("CrossValidationModels/SingleSp6predsCVtestEnvOnlySpatialEast",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(alpha,beta,draws,sumdraws,New,New.draws,logi,env.clim,env.clim.t)
    
    # 3. now only biotic predictors----------
    
    env.biot <- env[,c("MW2Prey","MW2CompR25")]
    env.biot.t <- env.t[,c("MW2Prey","MW2CompR25")]
    
    ## model
    n_env <- ncol(env.biot)
    n_sites <- nrow(env.biot)
    
    #greta
    alpha <- normal(0, 1)
    beta <- normal(0, 1, dim = n_env)
    
    # logit-linear model
    linear_predictor <- alpha + env.biot %*% beta
    p <- ilogit(linear_predictor)
    
    # distribution (likelihood) over observed values
    distribution(occupancy) <- bernoulli(p)
    
    #defining model
    m <- model(alpha,beta)
    
    #sampling 10.22
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 5000)
    
    sumdraws <- summary(draws)
    
    save(sumdraws,file=paste("CrossValidationModels/SingleSp6predsCVtrainBioticOnlySpatialEast",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    New <- alpha + env.biot.t %*% beta
    New.draws <- calculate(New,draws)
    
    logi <- summary(New.draws)
    logi$occur <- occupancy.t
    
    save(logi,file=paste("CrossValidationModels/SingleSp6predsCVtestBioticOnlySpatialEast",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(medi,dat,dat.test)
    
    
    
    rm(alpha,beta,draws,sumdraws,New,New.draws,logi,env,env.biot)
  }
}





runIndModPredsExamples <- function(sp.list,clim,bio1,bio2,type,metawebName){
  
  
  for (x in 1:length(sp.list)) {
    # for (x in 1:3) {
    print(paste("Start mod",sp.list[x],"at",Sys.time()))
    
    sp <- sp.list[x]
    
    load(paste("IndivModels/",sp,".siteIntersAllWebs.Rdata",sep=""))
    
    if (is.null(siteInters)) {print(paste("mod",sp.list[x],"failed to load"))}
    if (is.null(siteInters)) next
    
    t <- do.call("cbind",siteInters)
    
    datt <- merge(t,clim[,c(3,6:10)],by.x=0,by.y="PageName",all=F)
    colnames(datt)[23] <- 'foot'
    
    bio <- datt[,c(bio1,bio2)] 
    
    if (max(bio[,1])==0|max(bio[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio[,1])==0|max(bio[,2])==0) next
    
    cli <- datt[,c('b12','b6','b15','foot')]
    
    env <- cbind(bio,cli)
    env.sc <- apply(env,2,scale)
    
    #calc mean and sd used to scale original model
    mean.env <- apply(env,2,mean)
    sd.env <- apply(env,2,sd)
    occupancy <- datt[,"Occur"]
    
    #new data for responses to diff prey levels (25,50,75% prey)
    env.new <- data.frame(
      MW2Prey=rep(c(rep(quantile(env$MW2Prey)[2],100),
                    rep(quantile(env$MW2Prey)[3],100),
                    rep(quantile(env$MW2Prey)[4],100)),2),
      MW2CompR25=rep(0,600), 
      b12=rep(0,600),
      b6=c(rep(seq(min(env[,'b6']),max(env[,'b6']),length.out=100),3),rep(0,300)),
      b15=rep(0,600),
      foot=c(rep(0,300),rep(seq(min(env[,'foot']),max(env[,'foot']),length.out=100),3)))
    
    # scale new pred data using original scaling
    env.new.sc <- sweep(env.new,2,mean.env,FUN="-")
    env.new.sc <- sweep(env.new,2,sd.env,FUN="/")
    
    ## model
    n_env <- ncol(env)
    n_sites <- nrow(env)
    
    #greta
    alpha <- normal(0, 1)
    beta <- normal(0, 1, dim = n_env)
    
    linear_predictor <- alpha + env.sc %*% beta
    p <- ilogit(linear_predictor)
    
    distribution(occupancy) <- bernoulli(p)
    
    m <- model(alpha,beta)
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 5000)
    #predict - response curves
    New <- alpha + env.new.sc %*% beta
    New.draws <- calculate(New,draws)
    
    betas <- summary(draws)
    logi <- summary(New.draws)
    
    env.new$preds <- logi$statistics[,'Mean']
    env.new$q2p5 <- logi$quantiles[,'2.5%']
    env.new$q97p5 <- logi$quantiles[,'97.5%']
    
    save(betas,file=paste("ExampleModels/ExampleBetas",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    save(env.new,file=paste("ExampleModels/ExamplePredictions",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    
    rm(alpha,beta,draws,New,New.draws,logi)
  }
  
}





