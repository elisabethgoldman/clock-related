## read in percent methylation data for experimental dataset with necessary data imputed (see impute.R)
perc_meth_matrix <- readRDS("perc_meth_imputed.rds") 
## read in metadata file with chronological and otehr necessary info 
meta <- readDS("metadata.RDS")
chronologicalAge <- meta$Age

## read in two column tab-delimited file containing chrom_pos in one column and coef values in second column (note chrom_pos in column one should be same format as your rownames for your matrix of methylation datae.g., chr1_12345123)

sitesWithCoefs <- read.delim("mouseWLMTsitesCoefs.txt", sep="\t", header=TRUE)
## WLMT clock intercept
intercept <- 234.64 ## found at line 438 of Meer et al. Supplemental File 3 (Whole lifespan multi-tissue tab)

## subset percent methylation matrix (perc_imputed) to retain only sites in the clock that are also the perc meth matrix; format of rownames must match format of features ("chr_pos")
clockIDs <- sitesWithCoefs$Chrom_Position
clockValueInput <- perc_meth_imputed[rownames(perc_meth_imputed) %in% clockIDs,]
regWeights <- sitesWithCoefs$Weight

## calculate sum of weighted regression coefficients 
combine <- apply(clockValueInput, 2, function(x){x*regWeights)
predicted <- as.data.frame(colSums(combine))+intercept ## add intercept to totals

## examine strength of correlation with chorno age 
cor.test(predicted$`colSums(combine)`, chronologicalAge)
#calculate median absolute difference between predicted and observed 
mad <- median(abs(predicted-chronological Age))  #segregate by treatment group, then run pairwise one-sided two-sample t-tests, setting a priori expectation of which group (KO or experiemntal) is hypothesized to be "older" 
