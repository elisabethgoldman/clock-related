#!/usr/bin/env Rscript

library(bsseq)
library(boostme)
library(dplyr)
library(tidyselect)

methAll <- readRDS("unfilteredMethCounts.rds")
tcAll <- readRDS("unfilteredTCcounts.rds")


###convert data.frames to matrices

meth_impute <- as.matrix(methAll)
tc_impute <- as.matrix(tcAll)

##set NAs to zero (will be ignored by boostme)
tc_impute[is.na(tc_impute)] <- 0
meth_impute[is.na(meth_impute)] <- 0

##create objects chromosome and position to feed into bsseq 
cpg_windows <- rownames(meth_impute)
chr_position <- data.frame(do.call(rbind,strsplit(cpg_windows,"_",fixed = TRUE),quote=FALSE))

colnames(chr_position) <- c("chr", "start","end")
position <- as.character(chr_position$start)
num_position <- as.numeric(position)
chrom <- as.character(chr_position$chr)


sample_names <- colnames(meth_impute)
#sample_names2 <- gsub("\\./","",sample_names)

colnames(tc_impute) <- sample_names
colnames(meth_impute) <- sample_names
dim(tc_impute)
all.equal(rownames(meth_impute),rownames(tc_impute))## TRUE 
all.equal(colnames(meth_impute),colnames(tc_impute))## TRUE 

##create bsseq object

bs <- BSseq(M=data.matrix(meth_impute), Cov=data.matrix(tc_impute), sampleNames = sample_names, chr=chrom, pos=num_position)

##impute missing and low coverage values 

impute_sites <- boostme(bs,  randomCpGs = TRUE, minCov = 5, trainSize = 10000, testSize = 5000, validateSize = 5000, save = "boostMe_eval_metrics.txt", verbose = TRUE)

## Remove NaNs
perc_complete <- impute_sites[complete.cases(impute_sites),]
meth_complete <- meth_impute[rownames(meth_impute) %in% rownames(perc_complete),]
tc_complete <- tc_impute[rownames(tc_impute) %in% rownames(perc_complete),]

##save count data (unfiltered) except for having NaNs removed after imputation
saveRDS(meth_complete,file="completeImputedMethCounts.rds")
saveRDS(tc_complete,file="completeImputedTcCounts.rds")
saveRDS(perc_complete,file="completeImputedPercMethData.rds")


##Filter for >= 5X median coverage, between 10 and 90% median percent methylation (recommended but will need to see how much overlap you have with the clock site data) 
keep <- apply(tc_complete, 1, function(x){median(x,na.rm=TRUE)})>=5&apply(perc_complete,1,function(x)(median(x,na.rm=T)))>=0.1&apply(perc_complete, 1, function(x)(median(x,na.rm=T)))<=0.9
sum(keep)##total number of sites remaining after filtering
meth_filtered <- meth_complete[keep,]
tc_filtered <- tc_complete[keep,]
perc_filtered <- perc_complete[keep,]


#save filtered count data and perc meth matrix
saveRDS(meth_filtered, file="filteredCompleteImputedMethData.rds")
saveRDS(tc_filtered, file="filteredCompleteImputedTcData.rds")                                                                                                                              
saveRDS(perc_filtered,file="filteredCompleteImputedPercMethData.rds")
