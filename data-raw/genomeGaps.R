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
  gaps_hg19 <- keepSeqlevels(gaps_hg19, paste0("chr", c(1:22, "X", "Y")))
  save(gaps_hg19, file="~/Software/svpackages/svfilters/data/gaps_hg19.rda")
}
