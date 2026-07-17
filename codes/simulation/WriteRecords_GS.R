cat(" Writing records GS \n")
# NOTE: Training records are collected with a sliding-window process.
# Accumulation of the records starts in the burn-in year 'startTrainPop'.
# Once the burn-in period is over, the sliding-window process removes the oldest
# records.


if(year >= 11){
  
  ######------------------------------------------------
  ##>>>>>  pull the snps from male heterotic group
  ######------------------------------------------------
  
  if (year == 11) {
    train_yt1 = First_Clonal_Pheno
    train_yt1@fixEff <- as.integer(rep(paste0(year,1L),nInd(train_yt1)))
    
    train_yt2 = Second_Clonal_Pheno
    train_yt2@fixEff <- as.integer(rep(paste0(year,2L),nInd(train_yt2)))
    
    train_yt3 = Third_Clonal_Pheno
    train_yt3@fixEff <- as.integer(rep(paste0(year,3L),nInd(train_yt3)))
    
    
    trainPop = c(train_yt1, train_yt2, train_yt3)
    
  } else if (year >= 12 && year <= 20) {
    train_yt1 = First_Clonal_Pheno
    train_yt1@fixEff <- as.integer(rep(paste0(year,1L),nInd(train_yt1)))
    
    train_yt2 = Second_Clonal_Pheno
    train_yt2@fixEff <- as.integer(rep(paste0(year,2L),nInd(train_yt2)))
    
    train_yt3 = Third_Clonal_Pheno
    train_yt3@fixEff <- as.integer(rep(paste0(year,3L),nInd(train_yt3)))
    
    trainPop = c(trainPop, train_yt1, train_yt2, train_yt3)
    
  } else {
    train_yt1 = First_Clonal_Pheno
    train_yt1@fixEff <- as.integer(rep(paste0(year,1L),nInd(train_yt1)))
    
    train_yt2 = Second_Clonal_Pheno
    train_yt2@fixEff <- as.integer(rep(paste0(year,2L),nInd(train_yt2)))
    
    train_yt3 = Third_Clonal_Pheno
    train_yt3@fixEff <- as.integer(rep(paste0(year,3L),nInd(train_yt3)))

    #trainPop <- trainPop[-c((nInd(trainPop)+1-nInd(c(train_yt1, train_yt2, train_yt3))):nInd(trainPop))]
    trainPop <- trainPop[-c(1:nInd(c(train_yt1, train_yt2, train_yt3)))]
    trainPop = c(trainPop, train_yt1, train_yt2, train_yt3)
    
  }
  
}


