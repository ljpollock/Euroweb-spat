


runIndModPredsExamplesHypotheticalSpecies <- function(sp.list,clim,bio1,bio2,type,metawebName){
  
  
  for (i in 1:length(groups)) {
    gr <- groups[i]
    # mean parameters from species-level models
    params <- colMeans(meanparams.gold[grep(i,rownames(meanparams.gold)),])
    
    # find low, med, high prey 
    a <- s[grep(gr,s$ID),'Nprey']
    a.sc <- scale(a)
    
    cli <- clim[,c('b12','b6','b15','foot')]
    cli.sc <- apply(cli,2,scale)
    
    #new data for predictions (scaled)
    env.new <- data.frame(
      MW2Prey=rep(c(rep(quantile(a.sc)[2],100),
                    rep(quantile(a.sc)[3],100),
                    rep(quantile(a.sc)[4],100)),2),
      MW2CompR25=rep(0,600), 
      b12=rep(0,600),
      b6=c(rep(seq(min(cli.sc[,'b6']),max(cli.sc[,'b6']),length.out=100),3),rep(0,300)),
      b15=rep(0,600),
      foot=c(rep(0,300),rep(seq(min(cli.sc[,'foot']),max(cli.sc[,'foot']),length.out=100),3)))
    
    lin_pred <- params[1] + as.matrix(env.new) %*% params[c(2:7)] 
    pr <- exp(lin_pred) / (1 + exp(lin_pred))
    
    
    env.new$preds <- logi$statistics[,'Mean']
    env.new$q2p5 <- logi$quantiles[,'2.5%']
    env.new$q97p5 <- logi$quantiles[,'97.5%']
    
    save(betas,file=paste("ExampleModels/ExampleBetas",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    save(env.new,file=paste("ExampleModels/ExamplePredictions",metawebName,type,"SpSelect",sp,"Rdata",sep="."))
    
    
    rm(alpha,beta,draws,New,New.draws,logi)
  }
  
}




