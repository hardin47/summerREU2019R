
#  This function generates the size factors needed for normalization.

#  inputs
#     allCounts matrix
#     norm.samps are the columns to be normalized
#     norm.cond are the group levels


# norm.samps <- c("0.00_A", "0.00_B", "0.00_C", "0.35_A", "0.35_B", "0.35_C", 
#                "11.59_A", "11.59_B", "11.59_C", "20.40_A", "20.40_B", "20.40_C", 
#                "48.37_A", "48.37_B", "48.37_C", "129.96_A", "129.96_B", "129.96_C", 
#                "190.38_A", "190.38_B", "190.38_C")

# norm.cond <- as.factor(c("0.00", "0.00", "0.00", "0.35", "0.35", "0.35", "11.59", "11.59", "11.59", 
#                          "20.40",  "20.40",  "20.40", "48.37", "48.37", "48.37", "129.96", "129.96", 
#                          "129.96", "190.38", "190.38", "190.38"))

normfunc <- function(allCounts, norm.samps, norm.cond) {

# use DESeq2 to get size factors counts across all samples, all gene features
countsTable.all <- allCounts %>% dplyr::select(dplyr::one_of(norm.samps))

# define conditions

coldata.all <- as.data.frame(row.names =  colnames(countsTable.all), norm.cond)
dds.all <- DESeq2::DESeqDataSetFromMatrix(countData = countsTable.all, colData = coldata.all, design = ~norm.cond)
dds.all <- DESeq2::DESeq(dds.all)

# get results at the specified level of significance
res.all <- DESeq2::results(dds.all)

res.all <- as.data.frame(res.all) %>%
  dplyr::mutate(Geneid = allCounts$Geneid)

countsTable.all.normalized <- as.data.frame(DESeq2::counts(dds.all, normalized = TRUE)) %>%
  dplyr::mutate(Geneid = allCounts$Geneid, feature = allCounts$feature, genename = allCounts$genename)

# tidy countsTable normalized data
countsTable.all.normalized.tidy <-  countsTable.all.normalized %>%
  tidyr::gather(cond.samps, normCount, -Geneid, -feature, -genename) %>%
  dplyr::left_join(data.frame(cbind(norm.samps, norm.cond), stringsAsFactors = FALSE), 
                   by = c("cond.samps" = "norm.samps"))

list(normalized = countsTable.all.normalized, 
     normalizedT = countsTable.all.normalized.tidy, results = res.all)
}





