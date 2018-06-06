source("~/qiulab/cdg/cdg-functions.R")
x <- read.csv("~/qiulab/mskcc-sepsis-resistance/antibioticResistanceNorm.csv", header = TRUE, row.names = 1) #8 antibiotics table

x.snps <- read.csv("~/qiulab/mskcc-sepsis-resistance/processedSNPs.csv", header = TRUE, row.names = 1, stringsAsFactors = F) #snp table

yinhen.3 <- read.csv("~/qiulab/mskcc-sepsis-resistance/mutualReduced_3.csv", header = T, row.names = 1, na.strings = T)

yinhen.7 <- read.csv("~/qiulab/mskcc-sepsis-resistance/mutualReduced_7.csv", header = T, row.names = 1, na.strings = T)

x <- x[-which(rownames(x) == "PA7"),] #remove PA7 b/c it's another species

X <- list() #master list of all antibiotics
X.3 <- list() #master list of all SNPs reduced to 75% using MI algorithm
X.7 <- list() #master list of all SNPs reduced to 25% using MI algorithm

#make row names match tables x.snps and x
rownames(x.snps) <- sapply(strsplit(rownames(x.snps), "_"), "[", 1)
rownames(x.snps)[which(rownames(x.snps) == "PA01")] <- "PAO1"

x.new <- x.snps[row.names(x.snps) %in% row.names(x),] #reduced b/c database doesnt have all presented strains


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
  
  #Reduced list X.7
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


#loopOverList will cross-validate rows of each df in list & iterate 1000x
#rationale behind 1000x iteration is so that it will converge/stabalize & we capture the most frequently occuring highest ranked SNP
#note: some iterations in the same AB won't have same dim size b/c some small scoring SNPs don't show
ab_lists <- lapply(X.3, loopOverList) #list of all xgboost iterations of AB

ab_matrices <- lapply(ab_lists, function(x){ do.call(rbind, x)}) #ab_lists to 8 matrices

master <- do.call("rbind", ab_matrices)

ab_gain_mean <- aggregate(master$Gain, by = list(feature = master$Feature, AB = master$antibiotic_name), FUN = mean) #calculate mean by feature and ab name



##VISUALIZATION##
library(ggplot2)

#facet graph of the average of significant SNPs after 1000 iterations per AB
p <- ggplot(ab_gain_mean, aes(x = feature, y = x, label = feature)) + geom_point(aes(color = x)) + facet_wrap(~ AB) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=10), axis.text.x = element_text(angle=90, vjust=0.5, size=4))

p + geom_text(aes(label=ifelse(x > 0.2, as.character(feature),'')), hjust=0, vjust=0, size = 3) + labs(title="Mean of SNP Gains in 1000 runs per Antibiotic(8)")


#tree plotting

library(phytools)
cdg.tr <- read.tree(file = "~/qiulab/cdg/data/cdg-tree-mid.dnd")
cdg.names <- read.table(file = "~/qiulab/cdg/data/cdg.strains.txt3", row.names = 1, stringsAsFactors = F)
cdg.tr$tip.label <- cdg.names[cdg.tr$tip.label,1]

#make new x table re-ordered by phylo tree tips
#FIND MISMATCHING ROWS 
x.false <- which(rownames(x) %in% cdg.names[,1] == FALSE) #strain not in table x
cdg.false <- which(cdg.names[,1] %in% rownames(x)[-x.false] == FALSE) #strain not in table cdg.names
x.ab <- as.matrix(x[-x.false,])
x.ab <- rbind(x.ab, c(rep(NA,8))) #pad table x.ab w/ "W70332" to match cdg.names table
rownames(x.ab)[30] <- "W70332"
x.ab <- x.ab[match(cdg.tr$tip.label,row.names(x.ab)),] #reorder AB to tree tips

#reorder snps by tree tips
snp.reordered <- x.new[match(rownames(x.ab), rownames(x.new)),] 


#before using this function, please make sure snp.reordered variable exists
#this function is used to draw the phylogeny of 30 strains w/ causal SNPs alongside
displaySignificantSNPs <- function(AB){
  
  par(mar=c(5,1,4,0)) #fix margin to fit large graph
  ##DRAW PLOT SIZE & DISPLAY PHYLOGENY
  plot(cdg.tr, "p", use.edge.length = FALSE, x.lim = 70, cex = 0.8, main = toupper(AB)) #expand graph width w/ x.lim
  segments(rep(40, 30), 1:30, rep(40, 30) + x.ab[,AB], 1:30, lwd = 2, lty = "solid") #draw AB Resistance Index w/ segments function
  axis(1, at = c(38, 39, 40, 41,42,43), labels = c(-2,-1,0,1,2,3), cex.axis = 0.7)
  abline(v = 40)
  mtext("Antibiotic \nResistance \nIndex", at = 40, side = 1, line = 4, cex = 0.8)
  
  low <- ab_gain_mean[which(ab_gain_mean[,2] == AB & ab_gain_mean[,3] > 0.1 & ab_gain_mean[,3] < 0.2),1] #get snp names by AB & lt 0.2 but mt 0.1
  high <- ab_gain_mean[which(ab_gain_mean[,2] == AB & ab_gain_mean[,3] > 0.2),1] #get snp names by AB & mt 0.2

  name.ordered <- paste("SNP",sep="", sort(as.numeric(gsub(".*P","", c(low, high))))) #combine & sort low, high; this step is necessary for display
  
  mtext(name.ordered, side = 3, at = 50:(49+length(name.ordered)), cex = 0.6, las = 2, adj=0.2) #display SNP names
  
  for (i in 1:length(name.ordered)){
    label = name.ordered[i]
    col <- ifelse(match(label, high), "red", "black")
    text(49+i, 1:30, labels = snp.reordered[,label], cex = 0.7, col = col) #display SNPs lt 0.2
  }
}



#variable to change AB name
AB <- colnames(x.ab)[4] #can cycle through all 8
displaySignificantSNPs(AB)

#individual AB;
ggplot(ab_gain_mean[which(ab_gain_mean$AB == AB),], aes(x = feature, y = x, label = feature))+ geom_point(aes(color = x)) + geom_text(aes(label=ifelse(x > 0.15, as.character(feature),'')), hjust=0, vjust=0, size = 3) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=10), axis.text.x = element_text(angle=90, vjust=0.5, size=4), plot.margin = margin(4,7,4,10, "cm")) + labs(title=AB)


#create heatmap
#phylo.heatmap(cdg.tr, x)