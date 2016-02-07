library(svovarian)
library(svfilters)
data(lymph_ids)
dp <- projectOvarian()
pviews <- readRDS(file.path(dp["final_preprocess"], "pviews_hg19.rds"))
pviews <- pviews[, lymph_ids]
grl <- readRDS(file.path(dp["segment"], "grl_lymphoblasts_hg19.rds"))
g <- unlist(grl)
amps <- reduce(granges(g[g$seg.mean >= 1]))
dels <- reduce(granges(g[g$seg.mean <= -1]))
outs <- reduce(granges(germlineOutliers(pviews, NMAD=5)))
outs$type <- "outlier"
dels$type <- "deletion"
amps$type <- "amplicon"
lymphoblast_filters_hg19 <- GRangesList(list(amplicon=amps,
                                             deletion=dels,
                                             outlier=outs))
path <- "/dcl01/scharpf/data/svpackages/svfilters/data"
save(lymphoblast_filters_hg19, file=file.path(path, "lymphoblast_filters_hg19.rda"))


