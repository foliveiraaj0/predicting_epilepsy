setupEnvironment <- function() {
  localRepository <- "/home/fernandooliveira/Development/R/predicting_epilepsy"
  setwd(localRepository)
  source("http://bioconductor.org/biocLite.R")
  
  #biocLite("Biobase")
  #biocLite("GEOquery")
  #biocLite("genefilter")
  #install.packages("randomForest")
  #install.packages("e1071")
  #install.packages("DMwR")
  
  library(Biobase)
  library(GEOquery)
  library(genefilter)
  library(randomForest)
  library(e1071)
  library(DMwR)
  
  #GSE23703 = 506 samples
  # - GPL10806 = 483 samples
  # - GPL10807 = 23 samples
  
  #gsm <- getGEO(filename="GSE23703_family.soft", GSElimits = c(1, 3), GSEMatrix = FALSE)
}

summaryExpression <- function(expressionSet) {
  print("summarizing expression set ...")
  es <- exprs(expressionSet)
  #dim(es) -> 483 samples, each column is a sample. 180880 lines, each line is a Comparative Genomic Hybridization value 
  summary(as.vector(es))
}

resumeSampleStates <- function(expressionSet) {
  print("showing resume of sample states ...")
  pD <- phenoData(expressionSet)
  varMetadata(pD)
  table(expressionSet$characteristics_ch1.1)
}

loadRawData <- function() {
  print("loading raw data from remote server, this operation can take a while ...")
  gse <- getGEO("GSE23703", GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)
  expressionSet <- gse[[1]]
  return(expressionSet)
}

rowIQRs <- function(data) {
  rowQ(data,ceiling(0.75*ncol(data))) - rowQ(data,floor(0.25*ncol(data)))
}

plotData <- function(expressionSet) {
  #plot values without missing values or NA
  plot(rowMedians(expressionSet),rowIQRs(expressionSet),
       xlab='Median expression level',
       ylab='IQR expression level',
       main='Main Characteristics of Genes Expression Levels')
}

showData <- function(expressionSet) {
  hist(as.vector(expressionSet), breaks=80, prob=T, 
       xlab='Expression Levels', main='Histograma de Expression Levels')
  plotData(expressionSet)
}

filterMissingValues <- function(expressionSet) {
  print("filtering missing data and NA values, this opeartion can take a while ...")
  #searching for missing values or NA
  es <- exprs(expressionSet)
  valid <- complete.cases(es)
  expressionSet <- expressionSet[valid, ]
  return(expressionSet)
}

#not necessary for the tcc data
filterAffymetrixQualityControlProbeSets <- function(expressionSet) {
  print("filtering AFFimetrix quality control prob sets ...")
  idx <- grepl("AFFX", featureNames(expressionSet))
  expressionSet <- expressionSet[!idx,]
  return(expressionSet)
}

calculateIQRMean <- function(expressionSet) {
  print("calculating IQR mean, this opeartion can take a while ...")
  es <- exprs(expressionSet)
  iqrs <- c()
  for (feature in 1:nrow(es)) {
    iqrs[feature] <- IQR(es[feature,])
    # cat("iqrs[",feature,"] = ",iqrs[feature],"\n")
  }
  iqrsMean <- mean(iqrs)
  
  return(iqrsMean)
}

biocLoteFilterIQR <- function(expressionSet) {
  # not working!!
  print("applying biocLite IQR filter, this opeartion can take a while ...")
  es <- exprs(expressionSet)
  expressionSet <- nsFilter(expressionSet,
                   var.func=IQR,
                   var.cutoff=IQR(as.vector(es)/5),
                   feature.exclude="^AFFX")
  return(expressionSet)
}

filterIQR <- function(expressionSet) {
  print("applying IQR filter, this opeartion can take a while ...")
  es <- exprs(expressionSet)
  iqrsMean <- calculateIQRMean(expressionSet)
  idx <- c()
  for (feature in 1:nrow(es)) {
    idx[feature] <- IQR(es[feature,]) > iqrsMean
  }
  expressionSet <- expressionSet[idx,]
  return(expressionSet)
}

filterAnova <- function(expressionSet) {
  print("applying ANOVA filter, this opeartion can take a while ...")
  ff <- filterfun(Anova(expressionSet$characteristics_ch1.1, p = 0.0000001))
  es <- exprs(expressionSet)
  selGenes <- genefilter(es, ff)
  length(selGenes)
  expressionSet <- expressionSet[selGenes, ]
  return(expressionSet)
}

findMostImportantFeatures <- function(expressionSet) {
  print("finding the most important features, this operation can take a while ...")
  featureNames(expressionSet) <- make.names(featureNames(expressionSet))
  es <- exprs(expressionSet)
  dt <- data.frame(t(es), Mut = expressionSet$characteristics_ch1.1)
  rf <- randomForest(Mut ~ ., dt, importance = T)
  imp <- importance(rf)
  imp <- imp[, ncol(imp) - 1]
  rf.genes <- names(imp)[order(imp, decreasing = T)[1:30]]
  #sapply(rf.genes, function(g) tapply(dt[, g], dt$Mut, median))
  save(rf.genes, file = "mostImportantFeatures.Rdata")
  return(rf.genes)
}

filterData <- function(expressionSet) {
  print("starting filtering operation ...")
  expressionSet <- filterMissingValues(expressionSet)
  expressionSet <- filterIQR(expressionSet)
  expressionSet <- filterAffymetrixQualityControlProbeSets(expressionSet)
  expressionSet <- filterAnova(expressionSet)
  save(expressionSet, file = "expressionSet.Rdata")
  return(expressionSet)
}

run <- function() {
  setupEnvironment()
  if(file.exists("expressionSet.Rdata")) {
    print("filtered data found")
    load("expressionSet.Rdata")
    summaryExpression(expressionSet)
    load("mostImportantFeatures.Rdata")
  } else {
    expressionSet <- loadRawData()
    summaryExpression(expressionSet)
    expressionSet <- filterData(expressionSet)
    rf.genes <- findMostImportantFeatures(expressionSet)
  }
  
  genes <- lapply(rf.genes, function(x) substr(x, 2, nchar(x)))
  
  selectedFeatures <- featureNames(expressionSet) %in% genes
  
  expressionSet <- expressionSet[selectedFeatures]
  
  source("prediction.R")
  
  if(!file.exists("knn.Rdata") | !file.exists("svm.Rdata") | !file.exists("randomForest.Rdata")) {
    runModels(expressionSet)
  }
  
  compareResults('holdout')
  
}

compareResults <- function(method) {
  load("knn.Rdata")
  load("svm.Rdata")
  load("randomForest.Rdata")

  if(method == 'holdout') {
    
    models <- c(knn, svm, randomForest)
    
    for (modelIndex in 1:length(models)) {
      model <- models[[modelIndex]]
      iterations <- length(dimnames(model@foldResults)[[1]])
      verificationDataLength <- length(dimnames(model@foldResults)[[2]])
      learners <- length(dimnames(model@foldResults)[[3]])
      for (learnerIndex in 1:learners) {
        idx <- model@foldResults[1:iterations,1:verificationDataLength,learnerIndex:learnerIndex,1] == 100
        modelName <- dimnames(model@foldResults)[[3]][learnerIndex]
        accuracy <- sum(apply(idx, 1, function(x) sum(x))/verificationDataLength)/iterations
        cat(sprintf("%s = %.5f\n", modelName, accuracy))
        
      }
    }
  }
  else {
    all.trials <- join(knn, svm, randomForest, by = "variants")
    rankSystems(all.trials, top = 10, max = T)
  }

}

#lapply(b, function(x) substr(x, 2, nchar(x))) removing first character of all genes from b
#unlist(a, use.names=F) passing list to vector
