library(svfilters.hg18)
library(svcnvs)
germline_filters <- listGenomeFilters()
germline_filters <- germline_filters[-match("transcripts", names(germline_filters))]
save(germline_filters, file="svfilters.hg19/data/germline_filters.rda")
