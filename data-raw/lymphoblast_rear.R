library(svovarian)
library(svfilters)
library(svrearrange)
dp <- projectOvarian()
lymphoblast_rfiles <- file.path(dp["germline"], paste0("CGH", 1:10, "N.rds"))
lymphoblast_rear <- lapply(lymphoblast_rfiles, readRDS)
lb <- lapply(lymphoblast_rear, linkedBins)
lt <- lapply(lb, function(x) x$linked.to)
lb <- GRangesList(lb)
lt <- GRangesList(lt)
lb <- unlist(lb)
lt <- unlist(lt)
rear_intervals <- reduce(c(granges(lb), lt))
lymphoblast_rear_hg19 <- rear_intervals
save(lymphoblast_rear_hg19,
     file="/dcl01/scharpf/data/svpackages/svfilters/data/lymphoblast_rear_hg19.rda")





