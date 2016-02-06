library(svovarian)
library(svcnvs)
data(lymph_ids)
dp <- projectOvarian()
pviews <- readRDS(file.path(dp["final_preprocess"], "pviews_hg19.rds"))
pviews <- pviews[, lymph_ids]
grl <- readRDS(file.path(dp["segment"], "grl_lymphoblasts_hg19.rds"))
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
path <- "/dcl01/scharpf/data/svpackages/svfilters/data"
lymphoblast_filters_hg19 <- tracks
save(lymphoblast_filters_hg19, file=file.path(path, "lymphoblast_filters_hg19.rda"))


