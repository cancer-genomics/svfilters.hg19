library(normalblood)
library(svclasses)
library(svfilters)
dp <- projectBlood()
data(ids)
ids <- ids[-1] ##problems with first sample
rfiles <- file.path(dp["1somatic"], paste0(ids, ".rds"))
rglist <- lapply(rfiles, readRDS)
lb <- lapply(rglist, linkedBins)
lt <- lapply(lb, function(x) x$linked.to)
lb <- GRangesList(lb)
lt <- GRangesList(lt)
lb <- unlist(lb)
lt <- unlist(lt)
rear_intervals <- reduce(c(granges(lb), lt))
normalblood_rear_hg19 <- rear_intervals
save(normalblood_rear_hg19,
     file="/dcl01/scharpf/data/svpackages/svfilters/data/normalblood_rear_hg19.rda")
