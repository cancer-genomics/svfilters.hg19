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
cancer_genes <- read.csv(file.path("~/Software/svpackages/svfilters/inst/extdata",
                                   "drivers_2016-03-04.csv"),
                         header=TRUE, stringsAsFactors=FALSE)
drivers <- cancer_genes[, "Clinically.Significant"]
bio_sign <- cancer_genes[, "Biologically.Significant"]
is_driver <- drivers %in% bio_sign
remap <- read.delim("~/Downloads/grdoznk.txt", header=TRUE)
previous <- strsplit(remap$Previous.Symbols, ", ")
synonyms <- strsplit(remap$Synonyms, ", ")
## Some of the symbols are outdated
old_hgnc <- bio_sign [ ! bio_sign %in% tx$gene_name ]
old_hgnc2 <- gsub("ORF", "orf", old_hgnc)
old_hgnc2 <- gsub("SGK", "SgK", old_hgnc2)
new_hgnc <- rep(NA, length(old_hgnc2))
for(i in seq_along(old_hgnc)){
  cat(".")
  old <- old_hgnc2[i]
  ## Look first among previous symbols
  in_list <- sapply(previous, function(x, old) old %in% x, old=old)
  if(sum(in_list)==1){
    j <- which(in_list)
    new_hgnc[i] <- remap$Approved.Symbol[[j]]
    next()
  }
  if(substr(old, 1, 4) == "KIAA") next()
  ## Look under synonyms
  in_list <- sapply(synonyms, function(x, old) old %in% x, old=old)
  if(sum(in_list)==1){
    j <- which(in_list)
    new_hgnc[i] <- remap$Approved.Symbol[[j]]
    next()
  }
  if(old == "MLL4") {
    new_hgnc[i] <- "KMT2B"
    next()
  }
  if(old == "ZAK"){
    new_hgnc[i] <- "ZNF33A"
    next()
  }
  browser()
}
df <- data.frame(biol_sign=bio_sign)
index <- which(df$biol_sign %in% old_hgnc)
df$biol_sign[index] <- new_hgnc
df$is_clin_sign <- is_driver
df$previous_hgnc <- "NA"
df$previous_hgnc[index] <- old_hgnc2
write.csv(df, file="~/Software/svpackages/svfilters/inst/extdata/cancer_genes_2016-03-05.csv",
          row.names=FALSE)

