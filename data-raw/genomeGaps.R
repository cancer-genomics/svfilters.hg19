genomeGaps <- function(){
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome <- "hg19"
  mySession <- browserSession()
  genome(mySession) <- genome
  gaps <- getTable(ucscTableQuery(mySession, track="gap"))
  gaps.gr <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
                                         gaps$chromEnd),
                     type=gaps$type)
  gaps_hg19 <- gaps.gr
  ## no gaps available in mitochondria
  gaps_hg19 <- keepSeqlevels(gaps_hg19, paste0("chr", c(1:22, "X", "Y")), pruning.mode = "coarse") 
  seqinfo(gaps_hg19) <- seqinfo(Hsapiens)[seqlevels(gaps_hg19), ]
  gaps <- trim(gaps_hg19)
  save(gaps, file="/Users/dbruhm/Desktop/cancer-genomics/svfilters.hg19/data/gaps.rda")
}
library(BSgenome.Hsapiens.UCSC.hg18)
genome <- "hg18"
mySession <- browserSession()
genome(mySession) <- genome
gaps <- getTable(ucscTableQuery(mySession, track="gap"))
gaps.gr <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
                                       gaps$chromEnd),
                   type=gaps$type)
gaps_hg18 <- gaps.gr
## no gaps available in mitochondria
gaps_hg18 <- keepSeqlevels(gaps_hg18, paste0("chr", c(1:22, "X", "Y")))
seqinfo(gaps_hg18) <- seqinfo(Hsapiens)[seqlevels(gaps_hg18), ]
save(gaps_hg18, file="/dcl01/scharpf/data/svpackages/svfilters/data/gaps_hg18.rda")
