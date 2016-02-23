library(GenomicRanges)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg19)
slevels <- c(paste0("chr", 1:22), "chrX", "chrY")
bins <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg19)[slevels, ],
                   tilewidth=1000,
                   cut.last.tile.in.chrom=TRUE)
##
## Assembly gaps / GC extreme
##
gc <- GCcontent(Hsapiens, bins)
binAssemblyGaps <- bins[as.numeric(gc) < 0.1]
sum(width(binAssemblyGaps))/1e6 ## 234o Mb
binAssemblyGaps_hg19 <- binAssemblyGaps
save(binAssemblyGaps_hg19, file="../data/binAssemblyGaps_hg19.rda")

library(BSgenome.Hsapiens.UCSC.hg18)
slevels <- c(paste0("chr", 1:22), "chrX", "chrY")
bins <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg18)[slevels, ],
                   tilewidth=1000,
                   cut.last.tile.in.chrom=TRUE)
##
## Assembly gaps / GC extreme
##
gc <- GCcontent(Hsapiens, bins)
binAssemblyGaps_hg18 <- bins[as.numeric(gc) < 0.1]
sum(width(binAssemblyGaps_hg18))/1e6 ## 234o Mb
save(binAssemblyGaps_hg18, file="../data/binAssemblyGaps_hg18.rda")
