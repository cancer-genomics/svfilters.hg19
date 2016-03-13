library(svcnvs)
library(svfilters.hg19)
library(germline)
dp <- projectGermline()
grl <- readRDS(file.path(dp["segment"], "grl_hg19.rds"))
pviews <- readRDS(file.path(dp["final_preprocess"], "pviews_hg19.rds"))
g <- unlist(grl)
amps <- reduce(granges(g[g$seg.mean >= 1]))
dels <- reduce(granges(g[g$seg.mean <= -1]))
outs <- reduce(granges(germlineOutliers(pviews, NMAD=5)))
## get rid of chrY and chrM
dels <- reduce(keepSeqlevels(dels, paste0("chr", c(1:22, "X"))))
amps <- reduce(keepSeqlevels(amps, paste0("chr", c(1:22, "X"))))
outs <- reduce(keepSeqlevels(outs, paste0("chr", c(1:22, "X"))))
outs$type <- "outlier"
dels$type <- "deletion"
amps$type <- "amplicon"
## takes place of lymphoblast_filters_hg19
coverage_filters <- GRangesList(list(amplicon=amps, deletion=dels, outlier=outs))
save(coverage_filters,
     file="~/Software/svpackages/svfilters.hg19/data/coverage_filters.rda")
