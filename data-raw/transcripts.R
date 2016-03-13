library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
## Get all transcripts defined in Ensembl (version 75):
tx <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"))
library(GenomeInfoDb)
seqlevelsStyle(tx) <- "UCSC"
tx <- keepSeqlevels(tx, paste0("chr", c(1:22, "X", "Y", "M")))
tx <- sort(tx)
si <- seqinfo(tx)
genome(si) <- "hg19"
seqinfo(tx) <- si
df <- read.csv("~/Software/svpackages/svfilters/inst/extdata/cancer_genes_2016-03-05.csv",
               stringsAsFactors=FALSE)
df <- df[!is.na(df$biol_sign), ]
cancer.genes <- df$biol_sign[df$is_clin_sign]
tx$cancer_connection <- tx$gene_name %in% cancer.genes
tx$biol_sign <- tx$gene_name %in% df$biol_sign
## gene_name has the approved HGNC gene name
transcripts <- tx
save(transcripts, file="~/Software/svpackages/svfilters.hg19/data/transcripts.rda")
