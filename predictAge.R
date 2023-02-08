## read in percent methylation data for experimental dataset with necessary data imputed (see imputeMethylation.R)
perc_meth_imputed <- readRDS("filteredCompleteImputedPercMethData.rds") 
## read in metadata file with chronological and other necessary info 
meta <- readRDS("metadata.rds")
chronologicalAge <- meta$Age

## read in two column tab-delimited file containing chrom_pos in one column and coef values in second column (note chrom_pos in column one should be same format as your rownames for your matrix of methylation datae.g., chr1_12345123)

## For WLMT clock (Meer et al. 2018):
sitesWithCoefs <- read.delim("mouseWLMTsitesCoefs.txt", sep="\t", header=TRUE)
## WLMT clock intercept
intercept <- 234.64 ## found at line 438 of Meer et al. Supplemental File 3 (Whole lifespan multi-tissue tab)

## subset percent methylation matrix (perc_meth_imputed) to retain only sites part of the clock that are also in the perc meth matrix; format of rownames must match format across data obejcts or no  matches wil be made 

## Get loci for clock cpgs
clockIDs <- sitesWithCoefs$Chrom_Position
##subset 
clockInputs <- perc_meth_imputed[rownames(perc_meth_imputed) %in% clockIDs,]
regWeights <- sitesWithCoefs$Weight

## run weighted regression 
combine <- apply(clockInputs, 2, function(x){x*regWeights)
predicted <- as.data.frame(colSums(combine))+intercept ##sum to get age and add intercept == predicted age

## examine strength of correlation with chorno age 
cor.test(predicted$`colSums(combine)`, chronologicalAge)
#calculate median absolute difference between predicted and observed groups, use this difference statistic to run statistical tests 
mad <- median(abs(predicted - chronologicalAge))  #segregate data by treatment group ->  run pairwise, one-sided two-sample t-tests, setting a priori expectation of which group (KO or experiemntal) is hypothesized to be 'biologically' older  
