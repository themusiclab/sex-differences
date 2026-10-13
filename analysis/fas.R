### Exploratory and confirmatory factor analysis
require(corrplot)
require(psych)
require(qgraph)



e.fa <- function(d, loading.threshold, subsample=NULL){                   
  #corrplot(mixedCor(d,ncat=7), addCoef.col="grey",type="lower")
  if(!is.null(subsample)){
    d.tmp <- subsample
   }
  else{d.tmp <- d}
  cor.mat <- mixedCor(d.tmp,ncat = 7)$rho
  f1 <- fa(cor.mat,nfactors=1,cor="mixed",ncat=7)
  bad.vars <- NULL
  while(any(f1$loadings < loading.threshold)){
    bad.vars <- c(bad.vars, dimnames(f1$loadings)[[1]][f1$loadings < loading.threshold])
    print(bad.vars)
    d.tmp <- d.tmp[, !names(d.tmp) %in% bad.vars]
    cor.mat <- mixedCor(d.tmp,ncat = 7)$rho
    f1 <- fa(cor.mat,nfactors=1,cor="mixed", ncat=7)
    print(f1$loadings)
    }
  f1
 }

c.fa <- function(d, Bayesian=F, MCMC_sample=500, check_ordinal=T,missing="pairwise"){
  require(lavaan)
  m1 <- paste("f1 =~ ",paste(names(d)[!names(d) %in% "user_id"],collapse = "+"))
  if(check_ordinal==T){fit1 <- cfa(m1, data = d, missing=missing, ordered = is.ordinal(d))}
  else{
      fit1 <- cfa(m1, data = d, missing=missing,ordered = NULL)
  }
  if(Bayesian == T){
    require(blavaan)
    fit1 <- bcfa(m1, data = d, burnin = 2000,sample = MCMC_sample, n.chains = 2, target="stan",bcontrol=list(cores=4))
  }
  fit1
}

compute_fscores <- function(data=dat, indicator_vars, var_name, check_ordinal=T, missing="pairwise"){
    require(tidyverse)
    c.dat <- data[,c("user_id",indicator_vars)]
    rownames(c.dat) <- c.dat$user_id
    if(missing!="fiml"){
        #c.dat <- na.omit(c.dat)
        cf_model <- c.fa(c.dat, missing=missing,check_ordinal=check_ordinal)
    } else{
        cf_model <- c.fa(c.dat, missing = missing, check_ordinal=check_ordinal)
    }
    c.dat$fscores <- as.numeric(lavPredict(cf_model))
    colnames(c.dat)[colnames(c.dat)=="fscores"] <- var_name
    #c.dat <- rename(c.dat, var_name = fscores)
    c.dat$user_id <- rownames(c.dat)
    data <- merge(data, c.dat, all=T)
}
