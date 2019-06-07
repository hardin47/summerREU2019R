

#  This function takes DESeq2 results, uses the negative binomial to generate a different
#  count dataset, and outputs the simulated count & tidy datasets

# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory-behind-deseq2

parBS <- function(dds, allCounts, norm.samps, norm.cond){
  
  countsTable.all <- allCounts %>% dplyr::select(dplyr::one_of(norm.samps))
  
# define conditions
  
#  coldata.all <- as.data.frame(row.names =  colnames(countsTable.all), norm.cond)
#  dds <- DESeqDataSetFromMatrix(countData = countsTable.all, colData = coldata.all, design = ~norm.cond)
#  dds <- DESeq(dds)
  
  dds.muij <- SummarizedExperiment::assays(dds)[["mu"]]
  dds.size <- as.numeric(unlist(1/data.frame(disp = DESeq2::dispersions(dds))))

  ngenes <- dim(dds)[1]
  nsamp <- dim(dds)[2]
  newCounts <- matrix(nrow = ngenes, ncol=nsamp)
  for(i in 1:nsamp){
    mus <- unlist(dds.muij[,i])
     newCounts[,i] <- unlist(list(n=1, mu = mus, size = dds.size) %>% purrr::pmap(rnbinom))
  }
  colnames(newCounts) <- norm.samps
  
  newCounts <- allCounts %>%
    dplyr::select(Geneid, feature, Geneid1, Geneid2, feature1, feature2,
                  start.gene, end.gene, start.bnum, end.bnum, genename, bnum) %>%
    cbind(newCounts)


  # tidy countsTable normalized data
  newCounts.tidy <-  newCounts  %>%
    dplyr::select(-Geneid1, -Geneid2, -feature1, -feature2,
                  -start.gene, -end.gene, -start.bnum, -end.bnum, -bnum) %>%
    tidyr::gather(cond.samps, normCount, -Geneid, -feature, -genename) %>%
    dplyr::left_join(data.frame(cbind(norm.samps, norm.cond), stringsAsFactors = FALSE), 
                     by = c("cond.samps" = "norm.samps")) %>%
    dplyr::mutate(condition = as.numeric(norm.cond))
  
  list(counts = newCounts, countsT = newCounts.tidy)
}
