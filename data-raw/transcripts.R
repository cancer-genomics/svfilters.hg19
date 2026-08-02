##library(EnsDb.Hsapiens.v75)
##library(GenomeInfoDb)
##edb <- EnsDb.Hsapiens.v75
#### Get all transcripts defined in Ensembl (version 75):
##tx <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"))
##seqlevelsStyle(tx) <- "UCSC"
##tx <- keepSeqlevels(tx, paste0("chr", c(1:22, "X", "Y", "M")))
##tx <- sort(tx)
##si <- seqinfo(tx)
##genome(si) <- "hg19"
##seqinfo(tx) <- si
library(TxDb.Hsapiens.UCSC.hg19.refGene)
library(AnnotationHub)
tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.refGene)
ah <- AnnotationHub()
query(ah, "OrgDb")
hsapiens <- ah[["AH49582"]]
keytypes(hsapiens)
refseqid <- keys(hsapiens, "REFSEQ")
map <- select(hsapiens, refseqid, "SYMBOL", "REFSEQ")
table(tx$tx_name %in% map$REFSEQ)
## 66 of the transcripts are not in the map
tx$tx_name[! tx$tx_name %in% map$REFSEQ]
tx <- tx[ tx$tx_name %in% map$REFSEQ ]
symbols <- setNames(map$SYMBOL, map$REFSEQ)
symbols <- symbols[tx$tx_name]
tx$gene_name <- symbols

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
