

splineadjustment<-function(data,quantiles=c(0.25,0.5,0.75)){
  dm<-dimnames(data)
  dataq<-normalize.quantiles(data)
  dimnames(dataq)<-dm
  GeneMeans<-rowMeans(dataq,na.rm = TRUE)
  knots <- quantile(Means, p = quantiles)
  for(i in 1:ncol(dataq)){
    model <- lm (dataq[,i] ~ bs(GeneMeans, knots = knots),na.action=na.exclude)
    #model2<-smooth.spline(SangerMeans,SangerOverlap[,i],cv=TRUE)
    if(i==1){
      DataScreen<-dataq[,i]-fitted.values(model)+GeneMeans
    }else{
      DataScreen<-cbind(DataScreen,dataq[,i]-fitted.values(model)+GeneMeans)
    }
    
  }
  dimnames(DataScreen)<-dm
  return(DataScreen)
}


ComBatCP <- function (dat, batch, mod = NULL, par.prior = TRUE, 
                      mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam"),empBayes=TRUE) {
  ## make batch a factor and make a set of indicators for batch
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }
  if(length(dim(batch))>1){
    stop("This version of ComBat only allows one batch variable")
  }  ## to be updated soon!
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1+batch)  
  if (!is.null(ref.batch)){
    ## check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    cat("Using batch =",ref.batch, "as a reference batch (this batch won't change)\n")
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")
  
  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch  
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){
    mean.only=TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)
  
  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){
    check[ref] <- FALSE
  } ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  
  ## Number of covariates or covariate levels
  message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')
  
  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  
  ## Check for missing values
  NAs <- any(is.na(dat))
  if(NAs){
    message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
  ## print(dat[1:2,])
  
  ##Standardize Data across genes
  cat('Standardizing Data across genes\n')
  if (!NAs){
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else { 
    B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
  }
  
  ## change grand.mean for ref batch
  if(!is.null(ref.batch)){
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  }
  
  ## change var.pooled for ref batch
  if (!NAs){
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
    } else {
      var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
    }
  } else {
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)
    } else {
      var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
    }
  }
  
  stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
  }  
  s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
  
  ##Get regression batch effect parameters
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs){
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  } else{
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME
  }
  delta.hat <- NULL
  for (i in batches){
    if(mean.only==TRUE) {
      delta.hat <- rbind(delta.hat,rep(1,nrow(s.data))) 
    } else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
    }
  }
  
  if(empBayes){
  ##Find Priors
  gamma.bar <- rowMeans(gamma.hat)
  t2 <- rowVars(gamma.hat)
  a.prior <- apply(delta.hat, 1, sva:::aprior) # FIXME 
  b.prior <- apply(delta.hat, 1, sva:::bprior) # FIXME
  
  ## Plot empirical and parametric priors

  ## Find EB batch adjustments
  
  gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
  if (par.prior) {
    message("Finding parametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        gamma.star <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
        delta.star <- rep(1, nrow(s.data))
      }
      else {
        temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                             delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                             b.prior[i])
        gamma.star <- temp[1, ]
        delta.star <- temp[2, ]
      }
      list(gamma.star=gamma.star, delta.star=delta.star)
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  else {
    message("Finding nonparametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                         gamma.hat[i, ], delta.hat[i, ])
      list(gamma.star=temp[1,], delta.star=temp[2,])
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  }else{
    #no empirical bayes adjustment:
    gamma.star<-gamma.hat
    delta.star<-delta.hat
  }
  if(!is.null(ref.batch)){
    gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
    delta.star[ref,] <- 1  ## set reference batch variance equal to 1
  }
  
  ## Normalize the Data ###
  message("Adjusting the Data\n")
  
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j]))) # FIXME
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean # FIXME
  
  ## tiny change still exist when tested on bladder data
  ## total sum of change within each batch around 1e-15 
  ## (could be computational system error).  
  ## Do not change ref batch at all in reference version
  if(!is.null(ref.batch)){
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  
  return(list(correctedData=bayesdata,batchDesign=batch.design,gamma.star=gamma.star,delta.star=delta.star,varpool=var.pooled,stdmean=stand.mean))
}



BatchCorrection<-function(data1,data2,site1="Broad",site2="Sanger",CombatRes,stdPrior=TRUE,qcThresh=NULL,qcvalues1=NULL,qcvalues2=NULL){
  #need to make sure it's broad first then sanger for passing to adjustnewdata
  if(!is.null(qcThresh)){
    data1<-data1[,qcvalues1>=qcThresh]
    data2<-data2[,qcvalues2>=qcThresh]
  }
  site=c(rep(site1,ncol(data1)),rep(site2,ncol(data2)))
  adjusted<-AdjustNewData(data1,data2,CombatRes,site,stdPrior)
  return(adjusted)
}

AdjustNewData<-function(data1,data2,CombatRes,site,stdPrior=TRUE){
  blevels<-levels(as.factor(site))
  batch.design2<-CombatRes$batchDesign
  mean.star<-CombatRes$gamma.star
  var.star<-CombatRes$delta.star
  ngenes<-nrow(data1)
  mergedata<-cbind(data1,data2)
  sd1<-ncol(data1)
  
  dn<-dimnames(mergedata)
  mergedata<-normalize.quantiles(as.matrix(mergedata))
  dimnames(mergedata)<-dn
  data1<-mergedata[,colnames(data1)]
  data2<-mergedata[,colnames(data2)]
  usegenes<-intersect(rownames(data1),rownames(CombatRes$stdmean))
  data1<-data1[usegenes,]
  data2<-data2[usegenes,]
  ngenes<-length(usegenes)
  if(stdPrior){
    B_meanStd2<-matrix(CombatRes$stdmean[usegenes,1],nrow=ngenes,ncol=ncol(data1))
    S_meanStd2<-matrix(CombatRes$stdmean[usegenes,1],nrow=ngenes,ncol=ncol(data2))
    B_varStd2<-matrix(sqrt(CombatRes$varpool[usegenes,]),nrow=ngenes,ncol=ncol(data1))
    S_varStd2<-matrix(sqrt(CombatRes$varpool[usegenes,]),nrow=ngenes,ncol=ncol(data2))
  
    D1_all_std2<-(data1-B_meanStd2)/B_varStd2
    D2_all_std2<-(data2-S_meanStd2)/S_varStd2
  }else{
    grand.mean1<-rowMeans(data1)
    grand.mean2<-rowMeans(data2)
    stand.mean1 <- t(grand.mean1) %*% t(rep(1,ncol(data1)))
    stand.mean2 <- t(grand.mean2) %*% t(rep(1,ncol(data2)))
    sd1<-rowSds(data1)
    sd2<-rowSds(data2)
    B_varStd2<-t(sd1)%*%t(rep(1,ncol(data1)))
    S_varStd2<-t(sd2)%*%t(rep(1,ncol(data2)))
    
    D1_all_std2<-(data1-stand.mean1)/B_varStd2
    D2_all_std2<-(data2-stand.mean2)/S_varStd2
    
  }
  colnames(mean.star)<-rownames(CombatRes$stdmean)
  B_meanAll2<-matrix(mean.star[1,usegenes],nrow=ngenes,ncol=ncol(data1))
  S_meanAll2<-matrix(mean.star[2,usegenes],nrow=ngenes,ncol=ncol(data2))
  colnames(var.star)<-rownames(CombatRes$stdmean)
  B_var2<-sqrt(var.star[which(blevels=="Broad"),usegenes])%*%t(rep(1,ncol(data1)))
  S_var2<-sqrt(var.star[which(blevels=="Sanger"),usegenes])%*%t(rep(1,ncol(data2)))
  Broad_all_adjust2<-(D1_all_std2-B_meanAll2)/B_var2
  Sanger_all_adjust2<-(D2_all_std2-S_meanAll2)/S_var2
  
  Broad_all_adjust2<-(Broad_all_adjust2)*B_varStd2+B_meanStd2
  Sanger_all_adjust2<-(Sanger_all_adjust2)*S_varStd2+S_meanStd2
  #need to plot this data:
  alldata2<-cbind(Broad_all_adjust2,Sanger_all_adjust2)
  AdjData<-cbind(Broad_all_adjust2,Sanger_all_adjust2)
  dn<-dimnames(alldata2)
  alldata2<-normalize.quantiles(as.matrix(alldata2))
  dimnames(alldata2)<-dn
  return(list(qNorm=alldata2,NoNorm=AdjData))
  
}


RemovePC<-function(data,droppcanumber=1){
  if(sum(is.na(data))!=0){
    #Have NAs and need to impute missing values
    #data is genes x cell lines
    meanVals<-rowMeans(data,na.rm=TRUE)
    genesToimpute<-which(rowSums(is.na(data))!=0)
    for(i in 1:length(genesToimpute)){
      selcl<-which(is.na(data[genesToimpute[i],]))
      data[genesToimpute[i],selcl]<-meanVals[genesToimpute[i]]
    }
  }
  estpca<-prcomp(t(data),scale.=TRUE)
  npcas<-1:ncol(data)
  pcause<-npcas[!npcas%in%droppcanumber]
  df.denoised <- estpca$x[,pcause] %*% t(estpca$rotation[,pcause])
  df.denoised<-t(df.denoised)
  correctedData<-df.denoised*estpca$scale+estpca$center

    return(correctedData)

}


