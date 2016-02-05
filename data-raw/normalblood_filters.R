library(normalblood)
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
outs$seg.mean <- granges_copynumber(outs, pviews)
dels$seg.mean <- granges_copynumber(dels, pviews)
amps$seg.mean <- granges_copynumber(amps, pviews)
outs$type <- "outlier"
dels$type <- "deletion"
amps$type <- "amplicon"
tracks <- GRangesList(list(amplicon=amps, deletion=dels, outlier=outs))
path <- "~/Software/svpackages/svfilters/data"
normalblood_filters_hg19 <- tracks
save(normalblood_filters_hg19, file=file.path(path, "normalblood_filters_hg19.rda"))
