source("~/qiulab/cdg/cdg-functions.R")
x <- read.csv("~/qiulab/mskcc-sepsis-resistance/antibioticResistanceNorm.csv", header = TRUE, row.names = 1) #8 antibiotics table

x.snps <- read.csv("~/qiulab/mskcc-sepsis-resistance/processedSNPs.csv", header = TRUE, row.names = 1, stringsAsFactors = F) #snp table

yinhen.3 <- read.csv("~/qiulab/mskcc-sepsis-resistance/mutualReduced_3.csv", header = T, row.names = 1, na.strings = T)

yinhen.7 <- read.csv("~/qiulab/mskcc-sepsis-resistance/mutualReduced_7.csv", header = T, row.names = 1, na.strings = T)

x <- x[-which(rownames(x) == "PA7"),] #remove PA7 b/c it's another species

X <- list() #master list of all antibiotics
X.3 <- list() #master list of all SNPs reduced to 75%
X.7 <- list() #master list of all SNPs reduced to 25%

#process row names to match table x.snps and x
rownames(x.snps) <- sapply(strsplit(rownames(x.snps), "_"), "[", 1)
rownames(x.snps)[which(rownames(x.snps) == "PA01")] <- "PAO1"

x.new <- x.snps[row.names(x.snps) %in% row.names(x),]


for (i in 1:length(colnames(x))){
  #Not reduced list X
  antibiotic_name <- colnames(x)[i]
  data <- subset(x, select = c(antibiotic_name))
  merged <- merge(x.new, data, by = "row.names")
  X[[i]] <- as.matrix(merged[,c(126, 2:125)])
  rownames(X[[i]]) <- merged[,1]
  X[[i]] <- X[[i]][,-which(colSums(X[[i]], na.rm = TRUE) == 30)]
  
  #Reduced list X.3
  #Dataset is reduced using mutual information feature selection algorithm
  data_3 <- cbind(X[[i]][,1], yinhen.3)
  colnames(data_3)[1] <- colnames(X[[i]])[1]
  X.3[[i]] <- as.matrix(data_3)
  
  data_7 <- cbind(X[[i]][,1], yinhen.7)
  colnames(data_7)[1] <- colnames(X[[i]])[1]
  X.7[[i]] <- as.matrix(data_7)
}


head(runXg.cv(X[[1]]))
head(runXg.cv(X.3[[1]])) 
# head(runXg.cv(X.yinhen.7))

X.xg <- runXg.cv(X.yinhen.3)[,1:2]
X.xg.m <- melt(X.xg)
ggplot(X.tet25, aes(Gain, Feature)) + geom_tile(aes(fill=Gain), colour = "black")

n <- 100
loopOverX <- function(x){
  iter <- 0
  inner <- list()
  antibiotic_name <- colnames(x)[1]
  for (k in 1:n){
    data <- runXg.cv(x)[,1:2]
    iter <- iter + 1
    # iteration <- rep(iter, length(data$Feature))
    iteration <- iter
    rank <- seq(1, length(data$Feature))
    inner[[k]] <- cbind(data, iteration, rank, antibiotic_name)
  }
  return(inner)
}


master_lists <- lapply(X.3, loopOverX) #list of all xgboost iterations of antibiotics

master_matrix <- lapply(master_lists, function(x){ do.call(rbind, x)})

master <- do.call("rbind", master_matrix)

score <- rev(1:5)

for (i in 1:8){
  for (j in 1:100){
    cbind(score,master_matrix[[i]][which(master_matrix[[i]]$iteration == j)[1:5],])
  }
}
