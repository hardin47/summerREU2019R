

#  This file generates the size factors needed for normalization.

#  inputs
#     allCounts matrix
#     norm.samps are the columns to be normalized
#     norm.cond are the group levels
#     reference level for DE


# DE.samps <- c("0.00_A",  "0.00_C", "0.35_A",  "0.35_C", 
#                "190.38_A", "190.38_B", "190.38_C")

# DE.cond <- as.factor(c("0.00", "0.00", "0.00", "0.00", "190.38", "190.38", "190.38"))
 

# DEfunc(dp.countfilt, DE.samps, DE.cond, DE.ref = "0.00")

# DE.ref= "0.00"

DEfunc <- function(allCounts, DE.samps, DE.cond) {

DE.counts <- allCounts %>%
    dplyr::select(DE.samps)  
rownames(DE.counts) <- allCounts$Geneid

#genes that are DE across conditions
coldata.all <- as.data.frame(row.names =  colnames(DE.counts), DE.cond)
dds.all <- DESeq2::DESeqDataSetFromMatrix(countData = DE.counts, colData = coldata.all, design = ~DE.cond)
dds.all <- DESeq2::DESeq(dds.all)

dds.all

}

