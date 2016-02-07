library(normalblood)
library(svfilters)
library(svcnvs)
dp <- projectBlood()
data(ids)
ids <- ids[-1] ## problems with first sample
pviews <- readRDS(file.path(dp["final_preprocess"], "pviews_hg19.rds"))
pviews <- pviews[, ids]
paths(pviews) <- file.path(dp["3back"], basename(paths(pviews)))
grl <- lapply(list.files(dp["0cbs"], full.names=TRUE), readRDS)
grl <- GRangesList(grl)
g <- unlist(grl)
amps <- granges(g[g$seg.mean >= 1])
dels <- granges(g[g$seg.mean <= -1])
outs <- granges(germlineOutliers(pviews, NMAD=5))
## all are female
dels <- reduce(keepSeqlevels(dels, paste0("chr", c(1:22, "X"))))
amps <- reduce(keepSeqlevels(amps, paste0("chr", c(1:22, "X"))))
outs <- reduce(keepSeqlevels(outs, paste0("chr", c(1:22, "X"))))
outs$type <- "outlier"
dels$type <- "deletion"
amps$type <- "amplicon"
normalblood_filters_hg19 <- GRangesList(list(amplicon=amps, deletion=dels, outlier=outs))
path <- "/dcl01/scharpf/data/svpackages/svfilters/data"
save(normalblood_filters_hg19, file=file.path(path, "normalblood_filters_hg19.rda"))
