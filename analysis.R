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

n <- 1000
loopOverList <- function(x){
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


ab_lists <- lapply(X.3, loopOverList) #list of all xgboost iterations of antibiotics

ab_matrices <- lapply(ab_lists, function(x){ do.call(rbind, x)}) #list of 8 antibiotic matrices

master <- do.call("rbind", ab_matrices)

ab_gain_mean <- aggregate(master$Gain, by = list(feature = master$Feature, AB = master$antibiotic_name), FUN = mean)

p <- ggplot(ab_gain_mean, aes(x = feature, y = x, label = feature)) + geom_point(aes(color = x)) + facet_wrap(~ AB) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=10), axis.text.x = element_text(angle=90, vjust=0.5, size=4))

p + geom_text(aes(label=ifelse(x > 0.2, as.character(feature),'')), hjust=0, vjust=0, size = 3) + labs(title="Mean of SNP Gains in 1000 runs per Antibiotic(8)")
