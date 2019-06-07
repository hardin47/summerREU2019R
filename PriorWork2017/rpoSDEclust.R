
# This function runs the entire analysis from loading the data, normalizing
# DE analysis, and clustering.

fullDEclust <- function(allCounts, norm.samps, norm.cond, siglevel = 0.001, nclust = 6){

## Step 3. Normalize the full dataset & calculate p-values

dp.norm <- normfunc(allCounts, norm.samps, norm.cond)

## Step 5. Cluster the normalized counts that are significant for DE

# join normalized counts with significant DE
dp.normcount <- dp.norm$normalized %>% 
  dplyr::full_join(dp.norm$results, by = "Geneid") %>%
  dplyr::filter(padj < siglevel) %>%
  dplyr::select(-feature, -genename, -baseMean, -log2FoldChange, -lfcSE, -stat,
                -pvalue, -padj)
rownames(dp.normcount) <- dp.normcount[,"Geneid"]


dp.clustout <- clustfunc(dp.normcount, norm.samps, norm.cond, nclust = nclust)

dp.countclust <- dp.clustout$countclust
dp.medoids.k <- dp.clustout$medoids.k

list(countclust = dp.countclust, medoids.k = dp.medoids.k)

}
