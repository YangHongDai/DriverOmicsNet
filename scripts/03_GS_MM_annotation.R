# Load required libraries
library(UCSCXenaTools)
library(data.table)
library(tidyr)
library(dplyr)

# -----------------------------
# Step 1: Download Mutation Data
# -----------------------------
paraCohort <- 'TCGA Liver Cancer'
paraDatasets <- 'mc3_gene_level'

mut_TCGA <- XenaGenerate(subset = XenaHostNames == 'tcgaHub') %>%
  XenaFilter(filterCohorts = paraCohort) %>%
  XenaFilter(filterDatasets = paraDatasets)

XenaQuery(mut_TCGA) %>%
  XenaDownload(destdir = 'LIHC')

# -----------------------------
# Step 2: Process DNA Methylation Data
# -----------------------------
setwd('/Users/arthurdai/Desktop/Fighting/TCGA_integrative')

# Read annotation and ID mapping files
meth_anno <- fread('HM450.csv', data.table = FALSE)
meth_tcga <- read.csv('dna_meth_all_hg38', sep='\t')

# Load network genes (ensure dg is already defined)
gene_list <- dg
names(gene_list) <- 'Gene'

# Preprocess methylation ID mapping
tcga_sep <- separate_rows(meth_tcga, gene, sep = ",")
tcga_final <- tcga_sep[, c("X.id", "gene")]

# Filter based on network genes
meth_cancer <- subset(tcga_final, gene %in% gene_list$Gene)
meth_cancer_anno <- subset(meth_anno, Name %in% meth_cancer$X.id)

# Load and clean beta values
meth_cancer_beta <- fread('LIHC/meth_lihc.tsv', data.table = FALSE)
colnames(meth_cancer_beta) <- substr(gsub("-", ".", names(meth_cancer_beta)), 1, 15)
cancer_sample <- fread('LIHC/01_voomExpr_lihc.csv', data.table = FALSE)
meth_cancer_beta <- subset(meth_cancer_beta, select = colnames(meth_cancer_beta) %in% colnames(cancer_sample))

rownames(meth_cancer_beta) <- meth_cancer_beta[,1]
meth_cancer_beta[,1] <- NULL

# -----------------------------
# Step 3a: Promoter Region
# -----------------------------
meth_cancer_anno$TSS <- grepl("TSS", meth_cancer_anno$UCSC_RefGene_Group)
promoter_anno <- subset(meth_cancer_anno, TSS == TRUE)

promoter_beta <- subset(meth_cancer_beta, rownames(meth_cancer_beta) %in% promoter_anno$IlmnID)
promoter_beta$ID <- rownames(promoter_beta)
meth_cancer$ID <- meth_cancer$X.id

promoter_merged <- merge(promoter_beta, meth_cancer, by = "ID") %>%
  select(-ID, -X.id)

promoter_df <- promoter_merged %>%
  pivot_longer(cols = starts_with("TCGA"), names_to = "Sample", values_to = "Value") %>%
  group_by(gene, Sample) %>%
  summarise(Values = list(Value), .groups = 'drop') %>%
  pivot_wider(names_from = Sample, values_from = Values)

fwrite(promoter_df, file = "meth_lihc_matrix_promoter.csv")

# -----------------------------
# Step 3b: Non-CpG Island
# -----------------------------
meth_cancer_anno$CpG <- !grepl("Island", meth_cancer_anno$Relation_to_UCSC_CpG_Island)
island_anno <- subset(meth_cancer_anno, CpG == TRUE)

island_beta <- subset(meth_cancer_beta, rownames(meth_cancer_beta) %in% island_anno$IlmnID)
island_beta$ID <- rownames(island_beta)

island_merged <- merge(island_beta, meth_cancer, by = "ID") %>%
  select(-ID, -X.id)

island_df <- island_merged %>%
  pivot_longer(cols = starts_with("TCGA"), names_to = "Sample", values_to = "Value") %>%
  group_by(gene, Sample) %>%
  summarise(Values = list(Value), .groups = 'drop') %>%
  pivot_wider(names_from = Sample, values_from = Values)

fwrite(island_df, file = "meth_lihc_matrix_non_cpgisland.csv")

# -----------------------------
# Step 3c: Enhancer Region
# -----------------------------
meth_cancer_anno$Enhancer <- if_else(meth_cancer_anno$Enhancer == 'TRUE', 'TRUE', 'FALSE', missing = 'FALSE')
enhancer_anno <- subset(meth_cancer_anno, Enhancer == TRUE)

enhancer_beta <- subset(meth_cancer_beta, rownames(meth_cancer_beta) %in% enhancer_anno$IlmnID)
enhancer_beta$ID <- rownames(enhancer_beta)

enhancer_merged <- merge(enhancer_beta, meth_cancer, by = "ID") %>%
  select(-ID, -X.id)

enhancer_df <- enhancer_merged %>%
  pivot_longer(cols = starts_with("TCGA"), names_to = "Sample", values_to = "Value") %>%
  group_by(gene, Sample) %>%
  summarise(Values = list(Value), .groups = 'drop') %>%
  pivot_wider(names_from = Sample, values_from = Values)

fwrite(enhancer_df, file = "meth_lihc_matrix_enhancer.csv")