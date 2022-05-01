

#Prepares inputs for biotic interactions - selection based on species

runBioticInputsIndSpeciesWPrey <- function(MW1,MW2,MW4,site.by.species,sp.list) {
  
  for (x in 1:length(sp.list)) {
    
    sp <- sp.list[x] 
    
    
    siteInters <- vector("list", 17)
    names(siteInters) <- c("Occur",
                           "MW1Pred","MW1Prey","MW1CompR25","MW1ACompR25",
                           "MW2Pred","MW2Prey","MW2CompR25","MW2CompR90","MW2CompRA","MW2CompIntens","MW2CompRC","MW2ACompR25",
                           "MW4Pred","MW4Prey","MW4CompRA","MW4ACompRA")
    
    #choose sites
    O <- site.by.species[,sp]
    if (sum(O)<500) { 
      pres <- O[O==1] 
    } else {
      pres <- sample(O[O==1],500,replace=F) 
    }
    abs <- sample(O[O==0],1000,replace=F)
    locs <- c(pres,abs)
    siteInters[["Occur"]] <- locs
    
    
    
    #trim site by species to chosen sites
    Occur <- site.by.species[names(locs),]
    
    if (dim(Occur)[1]!= length(locs)) {print("warning:site by species does not match")}
    
    # metaweb 1
    PPmat <- MW1
    
    siteInters[["MW1Pred"]] <- (Occur %*% PPmat)[,sp]
    
    preypred <- t(PPmat)
    diag(preypred) <- 0
    
    siteInters[["MW1Prey"]] <- (Occur %*% preypred)[,sp]
    
    ShPrey <- PPmat %*% preypred
    diag(ShPrey) <- 0 
    
    #siteInters[["MW1CompIntens"]] <- Occur %*% ShPrey
    
    # Richness- only those species that share more than x% of prey (of total prey for species i)
    #divided out across rows.. row shares x% of its prey with column- Asymmetric
    cs <- ShPrey/Matrix::colSums(preypred)
    cs[is.na(cs) | cs<0.25] <- 0
    cs[cs>0.25] <- 1
    
    siteInters[["MW1CompR25"]] <- (Occur %*% t(cs))[,sp]
    rm(cs)
    
    #Apparent competition metrics--  
    ShPred <- preypred %*% PPmat 
    diag(ShPred) <- 0 
    
    cs <- ShPred/Matrix::colSums(PPmat)
    cs[is.na(cs) | cs<0.25] <- 0
    cs[cs>0.25] <- 1
    
    siteInters[["MW1ACompR25"]] <- (Occur %*% t(cs))[,sp]
    
    rm(PPmat,ShPred,cs,ShPrey)
    
    
    # metaweb 2
    PPmat <- MW2
    
    siteInters[["MW2Pred"]] <- (Occur %*% PPmat)[,sp]
    
    preypred <- t(PPmat)
    diag(preypred) <- 0
    
    siteInters[["MW2Prey"]] <- (Occur %*% preypred)[,sp]
    
    ShPrey <- PPmat %*% preypred
    diag(ShPrey) <- 0 
    
    siteInters[["MW2CompIntens"]] <- (Occur %*% ShPrey)[,sp]
    
    # Richness- only those species that share more than x% of prey (of total prey for species i)
    #divided out across rows.. row shares x% of its prey with column- Asymmetric
    cs <- ShPrey/Matrix::colSums(preypred)
    cs[is.na(cs) | cs<0.25] <- 0
    cs[cs>0.25] <- 1
    
    siteInters[["MW2CompR25"]] <- (Occur %*% t(cs))[,sp]
    rm(cs)
    
    cs <- ShPrey/Matrix::colSums(preypred)
    cs[is.na(cs) | cs<0.90] <- 0
    cs[cs>0.90] <- 1
    
    siteInters[["MW2CompR90"]] <- (Occur %*% t(cs))[,sp]
    rm(cs)
    
    ShPrey[ShPrey > 1] <- 1
    siteInters[["MW2CompRA"]] <- (Occur %*% ShPrey)[,sp]
    
    #filter again by co-occur
    tOccur <- t(Occur)
    ShSites <- tOccur %*% Occur
    rm(tOccur)
    ShSites[ShSites>1] <-1
    
    ShPreyR <- ShPrey * ShSites
    siteInters[["MW2CompRC"]] <- (Occur %*% ShPreyR)[,sp]
    
    #Apparent competition metrics--  
    ShPred <- preypred %*% PPmat 
    diag(ShPred) <- 0 
    
    cs <- ShPred/Matrix::colSums(PPmat)
    cs[is.na(cs) | cs<0.25] <- 0
    cs[cs>0.25] <- 1
    
    siteInters[["MW2ACompR25"]] <- (Occur %*% t(cs))[,sp]
    
    rm(PPmat,ShPred,cs,ShPrey)
    
    # metaweb 4
    PPmat <- MW4
    
    siteInters[["MW4Pred"]] <- (Occur %*% PPmat)[,sp]
    
    preypred <- t(PPmat)
    diag(preypred) <- 0
    
    siteInters[["MW4Prey"]] <- (Occur %*% preypred)[,sp]
    
    ShPrey <- PPmat %*% preypred
    diag(ShPrey) <- 0 
    
    ShPrey[ShPrey > 1] <- 1
    siteInters[["MW4CompRA"]] <- (Occur %*% ShPrey)[,sp]
    
    #Apparent competition metrics--  
    ShPred <- preypred %*% PPmat 
    diag(ShPred) <- 0 
    
    ShPred[ShPred > 1] <- 1
    
    siteInters[["MW4ACompRA"]] <- (Occur %*% ShPred)[,sp]
    
    save(siteInters,file=paste(sp,"siteIntersAllWebs.Rdata",sep="."))
    
    rm(siteInters)
  }
  
}


runBioticInputsIndSpeciesWPreyGBIFAmphs <- function(MW2,site.by.species,sp.list) {
  
  for (x in 1:length(sp.list)) {
    
    sp <- sp.list[x] 
    
    
    siteInters <- vector("list", 3)
    names(siteInters) <- c("Occur",
                           "MW2Prey","MW2CompR25")
    
    #choose sites
    O <- site.by.species[,sp]
    if (sum(O)<500) { 
      pres <- O[O==1] 
    } else {
      pres <- sample(O[O==1],500,replace=F) 
    }
    abs <- sample(O[O==0],1000,replace=F)
    locs <- c(pres,abs)
    siteInters[["Occur"]] <- locs
    
    
    
    #trim site by species to chosen sites
    Occur <- site.by.species[names(locs),]
    
    if (dim(Occur)[1]!= length(locs)) {print("warning:site by species does not match")}
    
    # metaweb 2
    PPmat <- MW2
    
    preypred <- t(PPmat)
    diag(preypred) <- 0
    
    siteInters[["MW2Prey"]] <- (Occur %*% preypred)[,sp]
    
    ShPrey <- PPmat %*% preypred
    diag(ShPrey) <- 0 
    
    # Richness- only those species that share more than x% of prey (of total prey for species i)
    #divided out across rows.. row shares x% of its prey with column- Asymmetric
    cs <- ShPrey/Matrix::colSums(preypred)
    cs[is.na(cs) | cs<0.25] <- 0
    cs[cs>0.25] <- 1
    
    siteInters[["MW2CompR25"]] <- (Occur %*% t(cs))[,sp]
    rm(cs)
    
    
    save(siteInters,file=paste(sp,"siteIntersMW2R25.Rdata",sep="."))
    
    rm(siteInters)
  }
  
}


# runs individual species models

runbioticPythoncomp <- function(bioclim,comptype,metawebName){
  
  species <- unique(bioclim$species) #11.02
  
  for (x in 1:length(species)) {
    #for (x in 6:length(species)) {  # testing
    print(paste("Start comp mod - R50",species[[x]],"at",Sys.time()))
    
    dat <- bioclim[bioclim$species==species[x],]
    
    env <- dat[,c("prey","b4","b12","b6","b15","foot")]
    
    if (comptype == "compRA") {
      comp <- dat$compRA
    }
    
    if (comptype == "compR50") {
      comp <- dat$compR50
    }
    
    if (comptype == "compRC") {
      comp <- dat$compRC
    }
    
    env <- cbind(env,comp)
    env <- apply(env,2,scale)
    
    occupancy <- dat[,c("Occur")]
    
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
    draws <- mcmc(m, warmup=1000,thin=10,n_samples = 10000)
    
    save(draws,file=paste("SingleSp7preds",metawebName,"R50","S1000",species[[x]],"Rdata",sep="."))
    
    rm(alpha,beta,env,occupancy,draws)   
  }
}



runIndModPython <- function(sp.list,clim,bio1,bio2,type,metawebName){
  
  
  for (x in 1:length(sp.list)) {
    # for (x in 1:3) {
    print(paste("Start mod",sp.list[x],"at",Sys.time()))
    
    sp <- sp.list[x]
    
    load(paste(sp,"siteIntersAllWebs.Rdata",sep="."))
    
    if (is.null(siteInters)) {print(paste("mod",sp.list[x],"failed to load"))}
    if (is.null(siteInters)) next
    
    t <- do.call("cbind",siteInters)
    
    dat <- merge(t,clim[,c(3,6:10)],by.x=0,by.y="PageName",all=F)
    colnames(dat)[23] <- 'foot'
    
    bio <- dat[,c(bio1,bio2)] 
    
    
    if (max(bio[,1])==0|max(bio[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio[,1])==0|max(bio[,2])==0) next
    ##
    cli <- dat[,c('b12','b6','b15','foot')]
    
    occupancy <- dat[,c("Occur")]
    
    env <- cbind(bio,cli)
    env <- apply(env,2,scale)
    
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
    
    save(sumdraws,file=paste("SingleSp6preds",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(alpha,beta,env,occupancy,draws,sumdraws)
  }
}






runIndModPythonLOGpreycomp <- function(sp.list,clim,bio1,bio2,type,metawebName){
  
  
  for (x in 1:length(sp.list)) {
    # for (x in 1:3) {
    print(paste("Start mod",sp.list[x],"at",Sys.time()))
    
    sp <- sp.list[x]
    
    load(paste(sp,"siteIntersAllWebs.Rdata",sep="."))
    
    if (is.null(siteInters)) {print(paste("mod",sp.list[x],"failed to load"))}
    if (is.null(siteInters)) next
    
    t <- do.call("cbind",siteInters)
    
    dat <- merge(t,clim[,c(3,6:10)],by.x=0,by.y="PageName",all=F)
    colnames(dat)[23] <- 'foot'
    
    bio <- dat[,c(bio1,bio2)] 
    
    if (max(bio[,1])==0|max(bio[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio[,1])==0|max(bio[,2])==0) next
    ##
    bio <- apply(bio,2,function(x) log(x+1))
    
    cli <- dat[,c('b12','b6','b15','foot')]
    
    occupancy <- dat[,c("Occur")]
    
    env <- cbind(bio,cli)
    env <- apply(env,2,scale)
    
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
    
    save(sumdraws,file=paste("SingleSp6predsLOG",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(alpha,beta,env,occupancy,draws,sumdraws)
  }
}





runIndModPythonLOGpreycompGBIF <- function(sp.list,clim,bio1,bio2,type,metawebName){
  
  
  for (x in 1:length(sp.list)) {
    # for (x in 1:3) {
    print(paste("Start mod",sp.list[x],"at",Sys.time()))
    
    sp <- sp.list[x]
    
    load(paste(sp,"siteIntersMW2R25.Rdata",sep="."))
    
    if (is.null(siteInters)) {print(paste("mod",sp.list[x],"failed to load"))}
    if (is.null(siteInters)) next
    
    t <- do.call("cbind",siteInters)
    
    dat <- merge(t,clim[,c(3,6:10)],by.x=0,by.y="PageName",all=F)
    colnames(dat)[9] <- 'foot'
    
    bio <- dat[,c(bio1,bio2)] 
    
    if (max(bio[,1])==0|max(bio[,2])==0) {print(paste("mod",sp.list[x],"has all zeros in biotic vars"))}
    if (max(bio[,1])==0|max(bio[,2])==0) next
    ##
    bio <- apply(bio,2,function(x) log(x+1))
    
    cli <- dat[,c('b12','b6','b15','foot')]
    
    occupancy <- dat[,c("Occur")]
    
    env <- cbind(bio,cli)
    env <- apply(env,2,scale)
    
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
    
    save(sumdraws,file=paste("SingleSp6predsLOGGBIF",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    rm(alpha,beta,env,occupancy,draws,sumdraws)
  }
}

