data(transcripts, package="svfilters.hg19")
drivers <- transcripts$gene_name[transcripts$cancer_connection]
drivers <- unique(drivers)
save(drivers, file="svfilters.hg19/data/drivers.rda")
