## functions for use with regularized regression models (ncvreg): -------

## 1st ncvreg function -------------------------------------------
# function to perform model selection with repeated CVs:
regulRegr_model_selection <- function(analysis.spec,mycluster){
  
  run.vector <- 1:analysis.spec$R
  
  analysis.name <- analysis.spec$analysis.name
  my.seed <- analysis.spec$myseed
  model.type <- analysis.spec$model.type
  
  ## ncvreg FULL model (clinical + molecular features): ------------------
  my.function.model.SCAD <- function(X,analysis.spec){
    require(ncvreg)
    analysis.name <- analysis.spec$analysis.name
    R <- analysis.spec$R
    K <- analysis.spec$K
    my.alpha <- analysis.spec$alpha
    my.weights <- analysis.spec$weights
    my.dfmax <- analysis.spec$dfmax
    dataset <- analysis.spec$dataset
    if(analysis.spec$scale.molecularData==TRUE){
      dataset$PredictorsMolecular <- as.data.frame(scale(data.matrix(dataset$PredictorsMolecular)))
    }
    my.seed <- analysis.spec$myseed
    set.seed(X+my.seed)
    model.type <- analysis.spec$model.type
    y.out = as.factor(dataset$Outcome)
    y = ifelse(y.out=="A",1,0)
    x = data.matrix(cbind(dataset$PredictorsClinical,dataset$PredictorsMolecular))
    
    # fitting SCAD model...................................................
    require(ncvreg)
    fit.success <- FALSE
    no.try <- 1
    while(!fit.success){
      modsel = try( cv.ncvreg(X = x,
                              y = y,
                              family = 'binomial',
                              penalty.factor = my.weights,
                              penalty = model.type,
                              max.iter=10000,
                              nfolds = K, 
                              alpha=my.alpha,
                              lambda.min=0.1,
                              dfmax=my.dfmax,
                              warn=TRUE),
                    silent=TRUE)
      fit.success <- ifelse(length(modsel)==1,FALSE,TRUE)
      if(!fit.success){
        print("trying again...")
        no.try <- no.try+1
        if(no.try > 10){ 
          print("fitting error!")
          no.try <- 0
          my.weights <- rep(1, length(my.weights))
          sink(file = paste(analysis.name, ".log"), append = TRUE)
          print(paste("Convergence problem in model fitting"))
          sink()
          #set.seed(format(Sys.time(), "%H%M%S"))
        }
      }
    } # end while !fit.success
    
    #a <- coef(modsel)
    #print(a[1:16])
    #sum(a[9:911])
    return(modsel)
    
  } # end my.function.model.SCAD
  
  ## ncvreg MOL model: ----------------------------
  my.function.model.MOL <- function(X,analysis.spec){
    require(ncvreg)
    analysis.name <- analysis.spec$analysis.name
    R <- analysis.spec$R
    K <- analysis.spec$K
    my.alpha <- analysis.spec$alpha
    my.weights <- analysis.spec$weights
    my.dfmax <- analysis.spec$dfmax
    dataset <- analysis.spec$dataset
    if(analysis.spec$scale.molecularData==TRUE){
      dataset$PredictorsMolecular <- as.data.frame(scale(data.matrix(dataset$PredictorsMolecular)))
    }
    my.seed <- analysis.spec$myseed
    set.seed(X+my.seed)
    model.type <- analysis.spec$model.type
    y.out = as.factor(dataset$Outcome)
    y = ifelse(y.out=="A",1,0)
    x = dataset$PredictorsMolecular
    
    # fitting MOL model........................................................
    require(ncvreg)
    fit.success <- FALSE
    no.try <- 1
    while(!fit.success){
      modsel = try( cv.ncvreg(X = x,
                              y = y,
                              family = 'binomial',
                              penalty.factor = rep(1,ncol(dataset$PredictorsMolecular)),
                              penalty = model.type,
                              max.iter=10000,
                              nfolds = K, 
                              alpha=my.alpha,
                              lambda.min=0.1,
                              dfmax=my.dfmax-ncol(dataset$PredictorsClinical),
                              warn=TRUE),
                    silent=TRUE)
      fit.success <- ifelse(length(modsel)==1,FALSE,TRUE)
      if(!fit.success){
        print("trying again...")
        no.try <- no.try+1
        if(no.try > 10){ 
          print("fitting error!")
          no.try <- 0
          set.seed(format(Sys.time(), "%H%M%S"))
        }
      }
    } # end while !fit.success
    
    #a <- coef(modsel)
    #print(a[1:16])
    #sum(a[9:911])
    return(modsel)
    
  } # end my.function.model.MOL
  
  ## run both models on cluster:
  clust.res.model.SCAD <- parLapply(cl=mycluster,X=run.vector,fun=my.function.model.SCAD,analysis.spec)
  clust.res.model.MOL  <- parLapply(cl=mycluster,X=run.vector,fun=my.function.model.MOL ,analysis.spec)
  
  ## save results:
  out <- list(clust.res.model.SCAD=clust.res.model.SCAD,
              clust.res.model.MOL=clust.res.model.MOL)
  # saveRDS(out,paste('../results/',analysis.name,'_modsel_',model.type,"_",my.seed,'.rds',sep=""))
  
  return(out)
  
} # end function regulRegr_model_selection

## function to evaluate model selection runs:-----------------------------------------------------------
evaluate_model_selections <- function(analysis.spec, ncvreg.fitted){
  
  analysis.name <- analysis.spec$analysis.name
  R <- analysis.spec$R
  K <- analysis.spec$K
  dataset <- analysis.spec$dataset
  if(analysis.spec$scale.molecularData==TRUE){
    dataset$PredictorsMolecular <- as.data.frame(scale(data.matrix(dataset$PredictorsMolecular)))
  }
  
  my.seed <- analysis.spec$myseed
  model.type <- analysis.spec$model.type
  pnames <- analysis.spec$proteinNames
  
  # modsel_results = readRDS(paste('../results/',analysis.name,'_modsel_',model.type,"_",my.seed,'.rds',sep=""))
  modsel.res.SCAD <- ncvreg.fitted$clust.res.model.SCAD
  modsel.res.MOL  <- ncvreg.fitted$clust.res.model.MOL
  
  model.fit.eval <- list()
  
  ##  analyse model fit for full model.SCAD: -------
 
  # -- Model Selection
  ## controls: cbind(coef(modsel_results[[1]]),modsel_results[[1]]$fit$beta[,modsel_results[[1]]$min])
  btas = do.call(cbind, lapply(modsel.res.SCAD, function(mod) as.numeric(coef(mod))))
  lams = do.call(c,     lapply(modsel.res.SCAD, function(mod) mod$lambda.min))
  
  ind = which(rowSums(btas)!=0)
  # no.clin.var <- ncol(dataset$PredictorsClinical)
  no.clin.var <- sum(ind[which(ind %in% ((1:ncol(dataset$PredictorsClinical))+1))]!=0)
  
  # checking if fit results can be evaluated:
  check <- "OK"
  
  # check the possibility that no proteins were selected:
  if(length(ind) == no.clin.var + 1) {
    print("No molecular predictors selected!")
    check <- "NOT OK: no molecular predictors"
  } 
  
  # check if intercept is part of the model (should always be):
  if(ind[1] != 1){
    print("CAUTION: no intercept model! Indexes might get wrong...")
    check <- "NOT OK: no intercept modelled"
  }
  
  if(check=="OK"){
    prot.signature = pnames[(ind[(no.clin.var+2):length(ind)]) - (ncol(dataset$PredictorsClinical)+1)] 
    
    beta.mean <- rowMeans(btas[ind[(no.clin.var+2):length(ind)],,drop=FALSE]) # to sort betas
    order.beta.mean <- order(beta.mean, decreasing = FALSE)
    
    # first is intercept, then clin.vars, then proteins
    par(mar=c(5,5,4,1)+.1)
    boxplot(t(btas[ind[(no.clin.var+2):length(ind)][order.beta.mean],,drop=FALSE]), # sort by incr beta
            names = prot.signature[order.beta.mean], 
            horizontal = TRUE, 
            las = 1,
            col = "khaki",
            main=paste("SCAD:",analysis.name))
    abline(v = 0, col = 'red', lwd = 1)
  
    # ---- Frequency:
    freq1u = rowMeans(btas != 0)
    freq.result <<- cbind(factor=c("int",
                                   colnames(dataset$PredictorsClinical)[((1:ncol(dataset$PredictorsClinical))+1) %in% ind], 
                                   prot.signature),
                          used=freq1u[freq1u!=0], # as.numeric?
                          beta=round(apply(btas[ind,],1,function(x){mean(x[x!=0])}),digits=7)) # as.numeric?
  
  
  
    model.fit.eval$prot.signature.SCAD=prot.signature
    model.fit.eval$freq.result.SCAD=freq.result
    model.fit.eval$all.betas.SCAD=btas[ind,] # 1st is intercept, then clin vars, then to length ind is protein signature
              
  } else {
    
    model.fit.eval$prot.signature.SCAD=NULL
    model.fit.eval$freq.result.SCAD=NULL
    model.fit.eval$all.betas.SCAD=NULL
            
  }
  
  ##  analyse model fit for model.MOL: ---------------------
  # -- Model Selection
  ## controls: cbind(coef(modsel_results[[1]]),modsel_results[[1]]$fit$beta[,modsel_results[[1]]$min])
  btas = do.call(cbind, lapply(modsel.res.MOL, function(mod) as.numeric(coef(mod))))
  lams = do.call(c,     lapply(modsel.res.MOL, function(mod) mod$lambda.min))
  
  ind = which(rowSums(btas)!=0)
  
  # ind.mol <- sort(ind, decreasing = TRUE)
  
  prot.signature = pnames[(ind[2:length(ind)]) - 1] 
 
  # checking if fit results can be evaluated:
  check <- "OK"
  
  # check the possibility that no proteins were selected:
  if(length(ind) <= no.clin.var + 1) {
    print("No molecular predictors selected!")
    check <- "NOT OK: no molecular predictors"
  } 
  
  # check if intercept is part of the model (should always be):
  if(ind[1] != 1){
    print("CAUTION: no intercept model! Indexes might get wrong...")
    check <- "NOT OK: no intercept modelled"
  }
  
  if(check=="OK"){
    
  beta.mean <- rowMeans(btas[ind[2:length(ind)],,drop=FALSE]) # to sort betas
  order.beta.mean <- order(beta.mean, decreasing = FALSE)
  
  # first is intercept, then proteins
  par(mar=c(5,5,4,1)+.1)
  boxplot(t(btas[ind[2:length(ind)][order.beta.mean],,drop=FALSE]), # sort by incr beta
          names = prot.signature[order.beta.mean],
          horizontal = TRUE, 
          las = 1,
          col = "khaki3",
          main=paste("MOL:",analysis.name))
  abline(v = 0, col = 'red', lwd = 1)
    
  # ---- Frequency:
  freq1u = rowMeans(btas != 0)
  freq.result <<- cbind(factor=c("int", prot.signature),
                        used=freq1u[freq1u!=0],
                        beta=round(apply(btas[ind,],1,function(x){mean(x[x!=0])}),digits=7))
                        #beta=btas[ind,1]) # 1st model beta of 500 repeats - 121021 ibv
                        
  model.fit.eval$prot.signature.MOL=prot.signature
  model.fit.eval$freq.result.MOL=freq.result
  model.fit.eval$all.betas.MOL=btas[ind,] # 1st is intercept, pos 2 to length ind is protein signature
  
  } else {
    
    model.fit.eval$prot.signature.MOL=NULL
    model.fit.eval$freq.result.MOL=NULL
    model.fit.eval$all.betas.MOL=NULL
    
  }
  
  
  ## output results for SCAD and MOL models:
  return(model.fit.eval)
  
} # end evaluate_model_selections function

## further aux functions: -------------------------------------------------------

## - Generates folds for binary classification, maintaining the ratio of cases to
##   controls in each fold.
bsplits = function(y, k) {
  
  caseids = which(y == 1)
  ctrlids = which(y == 0)
  
  shuffled_cases = sample(caseids)
  shuffled_ctrls = sample(ctrlids)
  
  case_splits = suppressWarnings(split(shuffled_cases, 1:k))
  ctrl_splits = suppressWarnings(split(shuffled_ctrls, 1:k))
  splits = mapply(c, case_splits, ctrl_splits)
  if (class(splits) == 'matrix') {splits = as.list(as.data.frame(splits))}
  
  foldid = rep(NA, length(y))
  for (i in 1:k) {
    foldid[splits[[i]]] = i
  } 
  
  return(foldid)
}

## function to deal with NAs in factors (which are not allowed for cv.nvreg):
replace.NAs <- function(myFactor){
  index.NA <- which(is.na(myFactor))
  if(is.factor(myFactor)){
    replacement <- sample(levels(myFactor),length(index.NA),replace=TRUE)
  } else {
    replacement <- rep(mean(myFactor,na.rm=TRUE),length(index.NA))
  }
  myFactor[index.NA] <- replacement
  return(myFactor)
}

## function to analyse performance (AUC, PPV, ...) from DCV performance runs:
analyse.perf <- function(perf.name,perf.values){
  #print(perf.name)
  ## analyse performance quality function:
  CI <-round(quantile(perf.values,probs=c(0.025,0.975),na.rm=TRUE),digits=3)
  meanCI <- c(round(mean(perf.values,na.rm=TRUE)-sd(perf.values,na.rm=TRUE)/sqrt(length(perf.values))*1.96,digits=3),
              round(mean(perf.values,na.rm=TRUE)+sd(perf.values,na.rm=TRUE)/sqrt(length(perf.values))*1.96,digits=3))
  #print(paste("mean:   ",round(mean(perf.values,na.rm=TRUE)*100,digits=2),"%",sep=""))
  #print(paste("CI:     ",CI[1],CI[2]))
  #print(paste("meanCI: ",meanCI[1],meanCI[2]))
  
  return(mean(perf.values,na.rm=TRUE))
}


# function for joining specific outcome and model
join.outcome.and.analysis <- function(out=NULL,ana=NULL){ # -> aux
  
  merged.data <- merge(x=data.frame(PID=out$PID,outc=out$outcome),                  # outcome: 288?
                       y=data.frame(PID=ana$PID,clin=ana$ClinData,mol=ana$MolData), # clin:    364 x 2, mol: 364 x 162?
                       sort=FALSE) 
  # merged:  263 x 166?
  colnames(merged.data)[(2+ncol(ana$ClinData)+1):ncol(merged.data)] <- colnames(ana$MolData)
  
  output.analysis <- list()
  output.analysis$outcome <- merged.data$outc
  output.analysis$ClinData <- merged.data[,3:(2+ncol(ana$ClinData)),drop=FALSE]
  output.analysis$MolData  <- merged.data[,(2+ncol(ana$ClinData)+1):ncol(merged.data)]
  
  return(output.analysis)
}

produce.analyses.to.run <- function(dataset.name,analysis.model,outcomes){
  
  my.analysis <- list() # to hold outcome and model, as well as parameters
  
  number.of.analyses <- length(analysis.model) * length(outcomes)
  
  print(paste("total number analyses if all have sufficient data:",number.of.analyses)) 
  # in ncvreg analysis: for each of these clin only, clin+proteins, proteins only
  # in RF analysis: only clin+proteins
  
  ana.counter <- 0
  my.analysis.names <- c()
  min.sample.size <- c()
  
  for(outcome.num in 1:length(outcomes)){
    
    for(analysis.num in 1:length(analysis.model)){
      
      ana.counter <- ana.counter + 1
      my.analysis.names[ana.counter] <- paste(outcomes[[outcome.num]]$name,"~",analysis.model[[analysis.num]]$name,sep="")
      
      my.analysis[[ana.counter]] <- join.outcome.and.analysis(out=outcomes[[outcome.num]],ana=analysis.model[[analysis.num]])
      
      my.analysis[[ana.counter]]$outcome.name <- paste(outcomes[[outcome.num]]$name)
      my.analysis[[ana.counter]]$model.name   <- paste(analysis.model[[analysis.num]]$name)
      
      #print(paste(ana.counter,my.analysis.names[ana.counter],":",table(my.analysis[[ana.counter]]$outcome)[1],table(my.analysis[[ana.counter]]$outcome)[2]))
      
      min.sample.size[ana.counter] <- min(c(table(my.analysis[[ana.counter]]$outcome)[1],table(my.analysis[[ana.counter]]$outcome)[2]))
      
    } # end analysis.model loop
    
  } # end outcomes loop
  
  ## select analyses with sufficient sample size:
  analyses.with.sufficient.samples <- which(min.sample.size >= 15) # use 13 for PedIBD_Uppsala
  print(paste(dataset.name,"models with sufficient numbers of samples (>15):",length(analyses.with.sufficient.samples)))
  
  analysis.data.to.run <- my.analysis[analyses.with.sufficient.samples]
  analysis.names.to.run <- my.analysis.names[analyses.with.sufficient.samples]
  
  return(list(dataset.name=dataset.name,
              analysis.data.to.run=analysis.data.to.run,
              analysis.names.to.run=analysis.names.to.run))
}

## function to cut long protein names in matrix
short.prot.names <- function(x) {
  a <- sub(".OID.*", "", colnames(x)) # remove OID 
  b <- gsub("^\\.+|\\.[^.]*$", "", a) # remove last part after dot (Uniprot ID)
  sh <- sub("mean.", "", b) # remove "mean" from averaged protein names
  # return(sh)
}

## function to cut long protein names in vector
short.prot.names.vec <- function(x) {
  a <- sub(".OID.*", "", x) # remove OID 
  b <- gsub("^\\.+|\\.[^.]*$", "", a) # remove last part after dot (Uniprot ID)
  sh <- sub("mean.", "", b) # remove "mean" from averaged protein names
  # return(sh)
}

# function to remove empty elements from list
list.remove.NULLs  <-  function(x) {
  x[unlist(lapply(x, length) != 0)]
}

## function to run ncvreg analysis -----------------------------------------------------------------------
run.ncvreg <- function(analysis.spec,mycluster,mycutoff) {
  
  print("Running ncvreg model fitting and variable selection...")
  # model fitting and variable selection:
  ncvreg.fit         <- regulRegr_model_selection(analysis.spec=analysis.syst,mycluster)
  ncvreg.fit.results <- evaluate_model_selections(analysis.spec=analysis.syst, 
                                                  ncvreg.fitted=ncvreg.fit)
  
  return(list(ncvreg.fit.res  = ncvreg.fit.results))
}