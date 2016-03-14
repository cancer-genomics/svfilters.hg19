## 
## Bins with low mappability already removed from 'bins1kb'
##
## - do we still need a low_mappability bins object?
##
mappability <- function(map, bins){
  ## summarize mappability by bin
  score <- rep(NA, length(bins))
  o <- findOverlaps(bins, map)
  if(length(o)==0) return(score)
  j <- subjectHits(o)
  i <- queryHits(o)
  subjectIndex <- split(j, i)
  w <- width(map)
  s <- map$score
  ## weight the mappability score by the width of the mappability interval
  avg <- sapply(subjectIndex, function(i, score, weight){
    (sum(score[i]*weight[i],na.rm=TRUE))/sum(weight[i], na.rm=TRUE)
  }, score=s, weight=w)
  score[as.integer(names(avg))] <- avg
  as.integer(score*1000)
}

library(rtracklayer)
library(svfilters.hg19)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biovizBase)
slevels <- c(paste0("chr", 1:22), "chrX", "chrY")
bins <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg19)[slevels, ],
                   tilewidth=1000,
                   cut.last.tile.in.chrom=TRUE)
path <- "/dcs01/rscharpf/seq/referenceGenome/hg19"
file <- file.path(path, "wgEncodeCrgMapabilityAlign100mer.bigWig")
map.bw <- import.bw(file)
low_mappability <-  map.bw[map.bw$score < 0.75]
low_mappability <- reduce(low_mappability, min.gapwidth=2e3)
pkgdir <- "/dcl01/scharpf/data/svpackages/svfilters.hg19/data"
save(low_mappability, file=file.path(pkgdir, "low_mappability.rda"))
##seqlevels(map.bw, force=TRUE) <- seqlevels(bins)
##map.gr <- subsetByOverlaps(map.bw, bins)
##map_score <- mappability(map.gr, bins)
######hist(map_score, breaks=1000)
##bins$map <- as.integer(map_score*1000)
##low_mappability <- bins[ bins$map < 750 & !is.na(bins$map) ] 
##pkgdir <- "/dcl01/scharpf/data/svpackages/svfilters.hg19/data"
##save(low_mappability, file=file.path(pkgdir, "low_mappability.rda"))
