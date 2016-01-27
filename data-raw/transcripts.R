getTranscripts <- function(){
  library(EnsDb.Hsapiens.v75)
  library(TxDb.Hsapiens.UCSC.hg19.refGene)
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
  cancer.genes <- read.csv(file.path("~/Software/svpackages/svfilters/inst/extdata",
                                     "master_4_3_2014.csv"),
                         skip=1, header=TRUE,
                         stringsAsFactors=FALSE)
  cancer.genes <- cancer.genes[["Gene.Symbol"]]
  cancer.genes[!cancer.genes %in% tx$gene_name]
  cancer.genes[match("MYCL1", cancer.genes)] <- "MYCL"
  tx$cancer_connection <- tx$gene_name %in% cancer.genes
  ## gene_name has the approved HGNC gene name
  tx_hg19 <- tx
  save(tx_hg19, file="~/Software/svpackages/svdata/data/tx_hg19.rda")
}
