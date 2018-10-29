
createVariants <- function() {
  
  vars <- list()
  
  #for complete results
  vars$randomForest <- list(ntree=c(500,750,100), mtry=c(5,15,30))
  vars$svm <- list(cost=c(1,100,500), gamma=c(0.01,0.001,0.0001))
  vars$knn <- list(k=c(3,5,7,11), norm=c(T,F))
  
  #for testing algorithm
  #vars$randomForest <- list(ntree=c(100), mtry=c(30))
  #vars$svm <- list(cost=c(500), gamma=c(0.0001))
  #vars$knn <- list(k=c(5,7,11), norm=c(T))
  
  return(vars)
}

genericModel <- function(form,train,test, learner, ...) {
    cat('=')
    tgt <- as.character(form[[2]])
    tgtCol <- which(colnames(train)==tgt)
    
    if (learner == 'knn') {
        pred <- kNN(form,
                      train,
                      test,
                      norm.stats=list(rowMedians(t(as.matrix(train[,-tgtCol]))),
                                      rowIQRs(t(as.matrix(train[,-tgtCol])))))
    }
    else {
       model <- do.call(learner,c(list(form,train),list(...)))
       pred <- if (learner != 'randomForest') predict(model,test) else predict(model,test,type='response')
    }
    
    c(accuracy=ifelse(pred == resp(form,test),100,0))
}

runModels <- function(expressionSet) {
  require(class,quietly=TRUE)
  require(randomForest,quietly=TRUE)
  require(e1071,quietly=TRUE)
  es <- exprs(expressionSet)
   # the dataset
  featureNames(expressionSet) <- make.names(featureNames(expressionSet))
  dt <- data.frame(t(exprs(expressionSet)),Mut=expressionSet$characteristics_ch1.1)
  DSs <- list(dataset(Mut ~ .,dt,'ALL'))
  vars = createVariants()
  # The learners to evaluate
  TODO <- c('knn','svm','randomForest')
  #TODO <- c('knn')
  for(td in TODO) {
    assign(td,
           experimentalComparison(
             DSs,
             c(
               do.call('variants',
                       c(list('genericModel',learner=td),
                         vars[[td]],
                         varsRootName=td))
             ),
             #loocvSettings(seed=1234,verbose=F)
             hldSettings(30,0.30,1234),
             itsInfo=T
           )
    )
    save(list=td,file=paste(td,'Rdata',sep='.'))
  }
}

kNN <- function(form, train, test, norm = T, norm.stats = NULL) {
  require(class, quietly = TRUE)
  tgtCol <- which(colnames(train) == as.character(form[[2]]))
  if (norm) {
    if (is.null(norm.stats))
      tmp <- scale(train[, -tgtCol], center = T, scale = T)
    else tmp <- scale(train[, -tgtCol], center = norm.stats[[1]],
                      scale = norm.stats[[2]])
    train[, -tgtCol] <- tmp
    ms <- attr(tmp, "scaled:center")
    ss <- attr(tmp, "scaled:scale")
    test[, -tgtCol] <- scale(test[, -tgtCol], center = ms,
                             scale = ss)
  }
  knn(train[, -tgtCol], test[, -tgtCol], train[, tgtCol])
}

testKNNFromClassLibrary <- function() {
  # library(class)
  # data(iris)
  # idx <- sample(1:nrow(iris), as.integer(0.7 * nrow(iris)))
  # tr <- iris[idx, ]
  # ts <- iris[-idx, ]
  # preds <- knn(tr[, -5], ts[, -5], tr[, 5], k = 3)
  # table(preds, ts[, 5])
  
  preds.norm <- kNN(Species ~ ., tr, ts, k = 3)
  table(preds.norm, ts[, 5])
}

# varsEnsembles <- function(tgt,train,test, varsSets, baseLearner,blPars, verb=F)  {
#     preds <- matrix(NA,ncol=length(varsSets),nrow=NROW(test))
#         
#     for(v in seq(along=varsSets)) {
#       
#         if (baseLearner=='knn') {
#             preds[,v] <- kNN(train[,-which(colnames(train)==tgt)], test[,-which(colnames(train)==tgt)],
#                            train[,tgt],blPars)
#         }
#         else {
#             m <- do.call(baseLearner, 
#                          c(list(as.formula(paste(tgt, paste(varsSets[[v]], collapse='+'),sep='~')),
#                                   train[,c(tgt,varsSets[[v]])]
#                                ),
#                           blPars)
#                         )
#                 
#             if (baseLearner == 'randomForest') {
#                 preds[,v] <- do.call('predict', list(m,test[,c(tgt,varsSets[[v]])], type='response'))
#             }
#             else {
#                 preds[,v] <- do.call('predict', list(m,test[,c(tgt,varsSets[[v]])]))
#             }
#         }
#     }
#     
#     ps <- apply(preds,1, function(x) levels(factor(x))[which.max(table(factor(x)))])
#     
#     ps <- factor(ps, levels=1:nlevels(train[,tgt]), labels=levels(train[,tgt]))
#     
#     if (verb) structure(ps,ensemblePreds=preds) else ps
# 
# }
