getSimulated <- function(snp1, snp2, nstrains, snps){
  n <- nstrains

  btwn <- snp1 * snp2
  
  covar.m <- rbind(c(1.0, snp1, snp2), c(snp1, 1.0, btwn), c(snp2, btwn, 1.0))
  
  sigmaEV <- eigen(covar.m)
  eps <- rnorm(n * ncol(sigmaEV$vectors))
  meps <- matrix(eps, ncol = n, byrow = TRUE)    
  meps <- sigmaEV$vectors %*% diag(sqrt(sigmaEV$values)) %*% meps #decomposition
  
  u <- t(pnorm(meps[2:nrow(meps),]))
  decretize.u <-qbinom(u,1,0.5)
  
  rand <- matrix(rnorm(nstrains*(snps-2)),nrow = nstrains)
  rand.u <- pnorm(rand)
  snp.rand <- qbinom(rand.u,1,0.5)
  x <- cbind(meps[1,], decretize.u, snp.rand)
  
  dimnames(x) <- list(c(), c("target", paste("feature", seq_len(ncol(x)-1), sep = ".")))
  
  return(x)
}

get1Simulated <- function(corco, nstrains, snps){
  r <- corco/10 # desired correlation coefficient
  sigma <- matrix(c(1,r,r,1), ncol=2) # var-covariance matrix
  s <- chol(sigma) #cholesky decomposition
  n <- nstrains 
  z <- s %*% matrix(rnorm(n*2), nrow=2)
  u <- pnorm(z[2,]) 
  snp.states <- qbinom(u, 1, 0.5)
  
  known <- t(rbind(z[1,],snp.states))
  rand <- pnorm(matrix(rnorm(nstrains*(snps-1)),nrow = nstrains))
  snp.rand<-qbinom(rand,1,0.5)
  x <- cbind(known, snp.rand)
  
  colnames(x)[1] <- "target"
  colnames(x)[-1] <- paste("feature", seq_len(ncol(x)-1), sep = ".")
  
  return(x)
}

#this cholesky method does not seem too stable compared to eigen decomposition but it is possible to generate multiple correlated variables using cholesky decomposition
getSimulated.chol <- function(snp1, snp2, nstrains, snps){
  n <- nstrains
  
  btwn <- snp1 * snp2
  
  covar.m <- matrix(c(1.0, snp1, snp2, snp1, 1.0, btwn, snp2, btwn, 1.0), nrow = 3) #cor-matrix 
  
  s <- chol(covar.m)
  z <- s %*% matrix(rnorm(n*3), nrow = 3)
  
  #decretize the 2 variables
  u <- t(pnorm(z[2:nrow(z),]))
  decretize.u <-qbinom(u,1,0.5)
  
  #generate some random discrete variables and combine w/ the 2 correlated variables
  rand <- matrix(rnorm(nstrains*(snps-2)),nrow = nstrains)
  rand.u <- pnorm(rand)
  snp.rand <- qbinom(rand.u,1,0.5)
  x <- cbind(z[1,], decretize.u, snp.rand)
  
  dimnames(x) <- list(c(), c("target", paste("feature", seq_len(ncol(x)-1), sep = ".")))
  
  return(x)
}

runXg <- function(x){
  require(xgboost)
  
  bst <- xgboost(data = x[,2:ncol(x)], label = x[,1], max.depth = 2, eta = .05, gamma = 0.3, nthread = 2, nround = 10, verbose = 0, eval_metric = "rmse")
  
  importance_matrix <- xgb.importance(model = bst, feature_names = colnames(x[,2:ncol(x)]))
  
  return(importance_matrix)
}

runXg.cv <- function(x){
  require(xgboost)

  n <- 1:nrow(x)
  p <- sample(n, nrow(x)*2/3) #partition
  
  train <- as.matrix(x[p,])
  test <- as.matrix(x[-p,])
  
  dtrain <- xgb.DMatrix(data = train[,2:ncol(train)], label = train[,1])
  dtest <- xgb.DMatrix(data = test[,2:ncol(test)], label = test[,1])
  
  watchlist <- list(train=dtrain, test=dtest)
  
  bst <- xgb.train(data=dtrain, nthread = 2, nround=10, watchlist=watchlist, eval.metric = "rmse", verbose = 0)
  
  importance_matrix <- xgb.importance(model = bst, feature_names = colnames(x[,2:ncol(x)]))
  
  return(importance_matrix)
  
} 


getPV <- function(df){

  p.value <- sapply(2:ncol(df), function(x){
    if((sum(df[,x], na.rm = T) > 1) & (!sum(df[,x], na.rm = T) == nrow(df)-1)){
      t.test(df[,1] ~ df[,x])$p.value
    }else if(sum(df[,x], na.rm = T) == nrow(df)-1){
      t.test(df[,1], mu = 1)$p.value
    }else {
      t.test(df[,1], mu= df[which(df[,x] == 1),][,1])$p.value
    }
  })

  p.value <- as.data.frame(p.value, row.names=colnames(df[,2:ncol(df)]))

  return(t(as.matrix(p.value)))
}


getTable <- function(x, importance_matrix, cor){
  core.f <- c("feature.1$", "feature.2$") #correlated features to target
  
  v.names <- c("feature", "SNPs", "strains", "cor_given", "cor_actual", "p.value", "gain", "cover", "rank")
  
  #use v.names to dynamically make numeric variables
  for (i in 1:length(v.names)){
    assign(v.names[i], numeric())
  }  
  
  n <- 0 #use to capture all length(x) of new_cor
  new_cor <- rep(cor, length(x)/length(cor))
  
  for (num in 1:length(x)) { 
    featNum <- NULL #this index needs to reset
    imp_matrix <- importance_matrix[[num]]
    x.iter <- x[[num]]
    x.pv <- getPV(x[[num]])
    
    #check presence of any core.f[i]
    #problem: if core.f[1] is not in imp_matrix but core.f[2] is, then it will be duplicated
    for (i in 1:2){ 
      if (length(grep(core.f[i], imp_matrix$Feature)) == 0){
        featNum <- c(featNum, i)
      } else {
        featNum <- c(featNum, grep(core.f[i], imp_matrix$Feature)) #get feature numbers{1,2} INDEX from each matrix  
      }
    }
    
    n <-  n + 1 #use to capture all length(x) of new_cor
    
    for(i in 1:2){
      f <- featNum[i]
      feature<-c(feature, imp_matrix$Feature[f])
      SNPs <-c(SNPs, 100)
      strains <- c(strains, nrow(x.iter))
      cor_given <- c(cor_given, new_cor[n])
      cor_actual <-c(cor_actual, cor(x.iter[,1],x.iter[,grep(paste(imp_matrix$Feature[f], "$", sep = ""), colnames(x.iter), perl = TRUE)]))
      p.value <- c(p.value, x.pv[,imp_matrix$Feature[f]])
      gain <- c(gain, imp_matrix$Gain[f])
      cover <- c(cover, imp_matrix$Cover[f])
      rank <- c(rank, f)
    }
  }
  
  x.master <- data.frame(feature, SNPs, strains, cor_given, cor_actual, p.value, gain, cover, rank)
  
  return(x.master)
}


