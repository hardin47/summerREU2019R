
#  This file generates two datasets: (1) counts, (2) tidy counts.  All the code is based on the very specific structure of the rpoS experiment with 24 samples.

#  inputs
#     filename


countfunc <- function(filename, norm.samps, norm.cond, countcutoff = 50){
  
allCounts <- read.csv(filename, header = T, sep = "\t", stringsAsFactors = FALSE)

bnum = "b[0-9]{4}"
allCounts$GeneidBackup = allCounts$Geneid
allCounts <- allCounts %>% tidyr::separate(GeneidBackup, c("feature", "rest"), sep="[:]")

allCounts <- allCounts %>% dplyr::select(Geneid, feature,                                           
                                  "0.00_A" = A0,
                                  "0.00_B" = B0,
                                  "0.00_C" = C0,
                                  "0.35_A" = A10..5,
                                  "0.35_B" = B10..5,
                                  "0.35_C" = C10..5,
                                  "11.59_A" = A5x10.5,
                                  "11.59_B" = B5x10.5,
                                  "11.59_C" = C5x10.5,
                                  "20.40_A" = A10..4,
                                  "20.40_B" = B10..4,
                                  "20.40_C" = C10..4,
                                  "48.37_A" = A10..3,
                                  "48.37_B" = B10..3,
                                  "48.37_C" = C10..3,
                                  "100.00_A" = A2537,
                                  "100.00_B" = B2537,
                                  "100.00_C" = C2537,
                                  "129.96_A" = A2x10..3,
                                  "129.96_B" = B2x10..3,
                                  "129.96_C" = C2x10..3,
                                  "190.38_A" = A5x10..3,
                                  "190.38_B" = B5x10..3,
                                  "190.38_C" = C5x10..3) %>%
  dplyr::select(Geneid, feature, dplyr::one_of(norm.samps))

# IGR's separate: 
# do start.bnum end.bnum start.genename end.genename
# left join igrs to allCounts
genename = ",[a-z]{3}[A-Z,]."
rna.name = ",rna[0-9].."
igr <- allCounts %>% dplyr::filter(feature %in% c("IGR", "AS_IGR"))
igr$GeneidBackup = igr$Geneid
igr <- igr %>% tidyr::separate(GeneidBackup, c("Geneid1", "Geneid2"), sep = "[/]")
igr$feature1 <- tidyr::separate(igr, Geneid1, c("feature1", "rest"), sep = "[,]")$feature1
igr$feature1 <- tidyr::separate(igr, feature1, c("rest", "feature1"), sep = "[()]")$feature1
igr$feature2 <- tidyr::separate(igr, Geneid2, c("feature2", "rest"), sep = "[,]")$feature2
igr$start.gene <- dplyr::case_when(
  igr$feature1 == "CDS" ~ stringr::str_extract(igr$Geneid1, genename),
  TRUE ~ stringr::str_extract(igr$Geneid1, rna.name))
igr$end.gene <- dplyr::case_when(
  igr$feature2 == "CDS" ~ stringr::str_extract(igr$Geneid2, genename),
  TRUE ~ stringr::str_extract(igr$Geneid2, rna.name))
igr$start.bnum <- dplyr::case_when(
  igr$feature1 == "CDS" ~ stringr::str_extract(igr$Geneid1, bnum),
  TRUE ~ "none")
igr$end.bnum <- dplyr::case_when(
  igr$feature2 == "CDS" ~ stringr::str_extract(igr$Geneid2, bnum),
  TRUE ~ "none")
igr <- igr %>% tidyr::separate(start.gene, into = c("comma", "start.gene"), sep = "[,]") %>% dplyr::select(-comma) %>% tidyr::separate(end.gene, into = c("comma", "end.gene"), sep = "[,]") %>% dplyr::select(-comma)
allCounts <- dplyr::full_join(igr, allCounts)

# CDS
# have bnum and genename columns
# left join to allCounts
genename = ":[a-z]{3}.."
cds <- allCounts %>% dplyr::filter(feature %in% c("AS_CDS", "CDS")) 
cds$genename <- stringr::str_extract(cds$Geneid, genename)
cds$bnum <- stringr::str_extract(cds$Geneid, bnum)
cds <- cds %>% tidyr::separate(genename, into = c("colon", "genename"), sep = ":") %>%
  dplyr::select(-colon)
allCounts <- dplyr::full_join(allCounts, cds)


#ncRNA
#ncRNA doesn't have bnums, but id's which we'll put in the genename column
rna.name = ":rna[0-9].."
rna <- allCounts %>% dplyr::filter(feature %in% c("ncRNA", "AS_ncRNA"))
rna$genename <- stringr::str_extract(rna$Geneid, rna.name)
rna <- rna %>% tidyr::separate(genename, into = c("colon", "genename"), sep = ":") %>%
  dplyr::select(-colon)
allCounts <- dplyr::full_join(allCounts, rna)

#rRNA
rRNA <- allCounts %>% dplyr::filter(feature %in% c("rRNA", "AS_rRNA"))
rRNA$genename <- stringr::str_extract(rRNA$Geneid, rna.name)
rRNA <- rRNA %>% tidyr::separate(genename, into = c("colon", "genename"), sep = ":") %>%
  dplyr::select(-colon)
allCounts <- dplyr::full_join(allCounts, rRNA)

#tRNA
tRNA <- allCounts %>% dplyr::filter(feature %in% c("tRNA", "AS_tRNA"))
tRNA$genename <- stringr::str_extract(tRNA$Geneid, rna.name)
tRNA <- tRNA %>% tidyr::separate(genename, into = c("colon", "genename"), sep = ":") %>%
  dplyr::select(-colon)
allCounts <- dplyr::full_join(tRNA, allCounts)

# remove the NA rows we just created by full_joining while adding the ncRNA, rRNA, tRNA genenames
allCounts <- dplyr::filter(allCounts, feature %in% c("IGR", "AS_IGR") | genename != "NA")

# make tidy data
countsTable.all.tidy <- allCounts %>%
  tidyr::gather(cond.samps, rawCount, -Geneid, -feature, -genename, -bnum, -Geneid1, -Geneid2, -feature1, -feature2, -start.gene, -end.gene, -start.bnum, - end.bnum)

# Keep only genes that have at least one sample with 50 counts
gene.big <- countsTable.all.tidy %>% dplyr::group_by(Geneid) %>% dplyr::filter(max(rawCount) >= countcutoff) %>% 
  dplyr::select(Geneid) %>% unique()

allCounts <- allCounts %>%
  dplyr::select(Geneid, feature, dplyr::one_of(norm.samps), Geneid1, Geneid2, feature1, feature2,
                start.gene, end.gene, start.bnum, end.bnum, genename, bnum)

countsTable.all.tidy <- countsTable.all.tidy %>%
  dplyr::filter(cond.samps %in% norm.samps)


allCounts <- allCounts %>% dplyr::filter(Geneid %in% gene.big$Geneid)
countsTable.all.tidy <- countsTable.all.tidy %>% dplyr::filter(Geneid %in% gene.big$Geneid)

list(counts = allCounts, countsT = countsTable.all.tidy)

}


