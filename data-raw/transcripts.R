
library(TxDb.Hsapiens.UCSC.hg19.refGene)
library(AnnotationHub)
tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.refGene)
ah <- AnnotationHub()
query(ah, "OrgDb")
hsapiens <- ah[["AH111575"]]
keytypes(hsapiens)
refseqid <- keys(hsapiens, "REFSEQ")
map <- select(hsapiens, refseqid, "SYMBOL", "REFSEQ")
table(tx$tx_name %in% map$REFSEQ)
## 433 of the transcripts are not in the map
tx$tx_name[! tx$tx_name %in% map$REFSEQ]
tx <- tx[ tx$tx_name %in% map$REFSEQ ]
symbols <- setNames(map$SYMBOL, map$REFSEQ)
symbols <- symbols[tx$tx_name]
tx$gene_name <- symbols

library(GenomeInfoDb)
seqlevelsStyle(tx) <- "UCSC"
tx <- keepSeqlevels(tx, paste0("chr", c(1:22, "X", "Y", "M")), pruning.mode = "coarse")
tx <- sort(tx)
si <- seqinfo(tx)
genome(si) <- "hg19"
seqinfo(tx) <- si

head(tx)
library(data.table)
library(dplyr)
library(readxl)

oncoKB_cancerGenes <- fread("data/oncokb_cancerGeneList.260106.tsv", 
                            header = T, sep = "\t", 
                            select = c("Hugo Symbol", "GRCh37 RefSeq","Gene Type", "Gene Aliases")) %>%
  rename(gene_name = "Hugo Symbol", tx_name = "GRCh37 RefSeq", 
         gene_type = "Gene Type", gene_alias = "Gene Aliases") %>%
  mutate(tx_name = gsub("\\..*","",tx_name)) %>%
  filter(grepl("ONCOGENE|TSG",gene_type))

#table(oncoKB_cancerGenes$gene_type)

oncoKB_cancerGenes_refSeq <- filter(oncoKB_cancerGenes, tx_name!="")


# Get list of genes involved in cancer (biol_sign) ------------------------
# Match genes based on refSeq ID
oncoKB_refseq <- oncoKB_cancerGenes_refSeq$tx_name
svfilter_refseq <- mcols(tx)$tx_name

refseq_match <- which(svfilter_refseq %in% oncoKB_refseq)


# If genes can't be matched based on refSeq id, match based on Hugo name.
# This may result in multiple matches for a single gene.
missing_refseq <- oncoKB_refseq[which(! oncoKB_refseq %in% svfilter_refseq)]
additional_Hugo <- filter(oncoKB_cancerGenes_refSeq, tx_name%in%missing_refseq) %>% 
  pull(gene_name)

oncoKB_Hugo <- filter(oncoKB_cancerGenes, tx_name=="" | tx_name%in%missing_refseq) %>%
  pull(gene_name)
svfilter_Hugo <- mcols(tx)$gene_name

hugo_match <- which(svfilter_Hugo %in% oncoKB_Hugo)

#tx$biol_sign <- FALSE
#tx$biol_sign[c(refseq_match, hugo_match)] <- TRUE

# Get list of clinically significant genes --------------------------------

# Summarize for each Gene, the row with the highest present Tx, Dx, Px, and R level (if present)
oncoKB_clinicallyRelevant <- 
  fread("data/oncokb_biomarker_Tx_Dx_Px.260112.txt", header = T) %>%
  select(Gene, Level) %>%
  filter(Gene!="Other Biomarkers") %>%
  mutate(Gene = gsub(" .*","",Gene)) %>%
  unique() %>%
  mutate(Level = if_else(Level %in% c("1", "2", "3","4"), paste0("Tx", Level), Level)) %>%
  mutate(
    clin_sig_class = sub("([A-Za-z]+).*", "\\1", Level),
    num      = as.integer(sub("[A-Za-z]+(\\d+)", "\\1", Level))
  ) %>%
  group_by(Gene,clin_sig_class) %>%
  filter(num == min(num)) %>%
  ungroup() %>%
  select(-clin_sig_class, -num) %>%
  arrange(Level) %>% 
  summarize(clinically_significant = paste0(Level,collapse = ","), .by = "Gene") %>%
  arrange(Gene)
  

# Convert GRanges to dataframe
tx_df <- as.data.frame(tx)

# Add column clinically_significant to tx
merged_df <- tx_df %>%
  left_join(oncoKB_clinicallyRelevant, by = join_by(gene_name == Gene))

tx2 <- makeGRangesFromDataFrame(merged_df, keep.extra.columns = TRUE)



# tx$cancer_connection <- FALSE
# tx$cancer_connection[clinicallyRelevant_match] <- TRUE


# Add list of cancer driver genes -----------------------------------------

# 825 genes in oncokb_cancerGeneList
oncokb_cancerGeneList <- pull(oncoKB_cancerGenes, gene_name) %>% unique()

# The 825 genes have 2,259 aliases, collectively
oncokb_gene_aliases <- pull(oncoKB_cancerGenes, gene_alias) %>%
  strsplit(split = ", ") %>% unlist() %>% unique()

oncokb_superset <- c(oncokb_cancerGeneList, oncokb_gene_aliases)


# PMID: 23539594
vogelstein2013_Science <- read_xlsx("data/vogelstein2013_Science.xlsx", sheet = 6, skip = 1, n_max = 125,   col_names = T) %>%
  pull(`Gene Symbol`) %>% unique()

# PMID: 24132290
kandoth2013_Nature <- read_xlsx("data/kandoth2013_Nature_TableS4.xlsx", sheet = 1, n_max = 128, col_names = T) %>%
  pull(Gene) %>% unique()

# PMID: 24390350
lawrence2014_Nature <- read_xlsx("data/lawrence2014_Nature_TableS2.xlsx", sheet = 1, skip = 1, col_names = T) %>%
  pull(gene) %>% unique()

# PMID: 29056346
martincorena2017_Cell <- read_xlsx("data/martincorena2017_Cell.xlsx", sheet = 1, col_names = T) %>%
  pull(Gene) %>% unique()

# PMID: 29625053
bailey2018_Cell <- read_xlsx("data/bailey2018_Cell.xlsx", sheet = 2, skip = 3, col_names = T) %>%
  pull(Gene) %>% unique()

# PMID: 32015527
dietlein2020_NatGen <- read_xlsx("data/dietlein2020_NatGen.xlsx", sheet = 5, col_names = T) %>%
  filter(`Literature Support` %in% c("level A","level B")) %>%
  filter(`Min FDR` < 0.05) %>%
  pull(Gene) %>% unique()

driver_genes <- c(vogelstein2013_Science,kandoth2013_Nature,lawrence2014_Nature,martincorena2017_Cell,bailey2018_Cell,dietlein2020_NatGen)

# union(oncokb_cancerGeneList, driver_genes) %>% length()
# intersect(oncokb_gene_aliases, driver_genes)


# missing_from_oncoKB <- setdiff(driver_genes, oncokb_superset)
# length(missing_from_oncoKB)

cancer_genes_all <- unique(c(driver_genes, oncokb_superset, oncoKB_clinicallyRelevant$Gene))
cancer_genes_match <- which(mcols(tx)$gene_name %in% cancer_genes_all)

# In total, there are 1,062 unique "cancer genes" in the tx object
length(unique(tx$gene_name[cancer_genes_match]))

tx2$cancer_gene <- FALSE
tx2$cancer_gene[cancer_genes_match] <- TRUE

# Save transcripts object -------------------------------------------------

transcripts <- tx2
save(transcripts, file="data/transcripts.rda")
