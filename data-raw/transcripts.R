getTranscripts <- function(){
  library(EnsDb.Hsapiens.v75)
  library(TxDb.Hsapiens.UCSC.hg19.refGene)
  library(bioMaRt)
  edb <- EnsDb.Hsapiens.v75
  ## Get all transcripts defined in Ensembl (version 75):
  tx <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"))
  library(GenomeInfoDb)
  seqlevelsStyle(tx) <- "UCSC"
  tx <- keepSeqlevels(tx, paste0("chr", c(1:22, "X", "Y", "M")))
  tx <- sort(tx)
  ## gene_name has the approved HGNC gene name
  tx_hg19 <- tx
  save(tx_hg19, file="~/Software/svpackages/svdata/data/tx_hg19.rda")
}
