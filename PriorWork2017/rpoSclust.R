

#  This file generates the size factors needed for normalization.

#  inputs
#     normalized, filtered data
#     number of clusters

# clustfunc(dp.normcount, nclust = 6)

clustfunc <- function(normcount, norm.samps, norm.cond, nclust) {
  
   
  #generate dissmilarity matrix based on spearman correlation
  dissS <- 1 - cor(t(subset(normcount, select = -Geneid)), method = "spearman")
  
  pam.output <- cluster::pam(dissS, nclust)
  thisCluster <- data.frame(cluster = pam.output$clustering)
  pam.medoids <- data.frame(medoids = rownames(pam.output$medoids))
  
  pam.clusters <- thisCluster %>%
    dplyr::mutate(Geneid = normcount$Geneid)
  
  # make tidy data and RpoS level
  normcountT <- normcount %>%
    tidyr::gather(cond.samps, normCount, -Geneid) %>%
    dplyr::left_join(data.frame(cbind(norm.samps, norm.cond), stringsAsFactors = FALSE), 
                     by = c("cond.samps" = "norm.samps")) %>%
    dplyr::mutate(RpoS = as.numeric(norm.cond)) %>%
    dplyr::group_by(Geneid, RpoS) %>%
    dplyr::mutate(levelMean = mean(normCount)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Geneid) %>%
    dplyr::mutate(scaler = max(levelMean[1], levelMean[length(table(cond.samps))]),  #first and last mean
           normCountScaled = normCount/scaler,
           meanScaled = levelMean/scaler) %>%
    dplyr::ungroup() 
  
  countclust <- dplyr::left_join(pam.clusters, normcountT, by = "Geneid") 
  medoids.k <- dplyr::filter(countclust, Geneid %in% pam.medoids$medoids)
  
  list(countclust = countclust, medoids.k = medoids.k)
  
}


  
  