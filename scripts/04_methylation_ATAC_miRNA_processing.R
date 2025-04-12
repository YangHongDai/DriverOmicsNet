### Load necessary libraries
library(UCSCXenaTools)
library(data.table)
library(dplyr)
library(tidyr)
library(ceRNAnetsim)
library(miRBaseConverter)
library(foreach)

### Download and process ATAC-seq data (promoter and enhancer regions)
paraCohort <- 'GDC Pan-Cancer'
paraDatasets <- 'TCGA_ATAC_peak_Log2Counts_dedup_promoter'

ATAC_TCGA <- XenaGenerate(subset = XenaHostNames == 'atacseqHub') %>%
  XenaFilter(filterCohorts = paraCohort) %>%
  XenaFilter(filterDatasets = paraDatasets)

XenaQuery(ATAC_TCGA) %>% XenaDownload(destdir = 'COAD/ATAC')

pan_atac <- fread('ATAC/TCGA_ATAC', data.table = FALSE)
rownames(pan_atac) <- pan_atac$sample
names(pan_atac) <- substr(gsub("\\-", ".", names(pan_atac)), 1, 15)
cancer_sample <- fread('LIHC/01_voomExpr_lihc.csv', data.table = FALSE)
cancer_sample.2 <- cancer_sample[, grepl('TCGA', colnames(cancer_sample))]
atac_df <- pan_atac[, colnames(pan_atac) %in% colnames(cancer_sample.2)]
atac_df$id <- rownames(atac_df)

# Process promoter ATAC-seq
gene_list <- dg
promoter_map <- fread('ATAC/Promoter_map', data.table = FALSE) %>% separate_rows(gene, sep = ",")
atac_promoter <- left_join(atac_df, promoter_map, by = 'id') %>%
  filter(gene %in% gene_list$Gene) %>%
  select(-one_of(names(promoter_map)[c(1, 3, 4, 5, 6)]))

df_atac_promoter <- atac_promoter %>%
  pivot_longer(cols = starts_with("TCGA"), names_to = "Sample", values_to = "Value") %>%
  group_by(gene, Sample) %>%
  summarise(Values = list(Value)) %>%
  pivot_wider(names_from = Sample, values_from = Values)

fwrite(df_atac_promoter, file = "atac_lihc_promoter.csv")

# Process enhancer ATAC-seq
enhancer_map <- fread('ATAC/Enhancer_map', data.table = FALSE) %>% separate_rows(gene, sep = ",")
atac_enhancer <- left_join(atac_df, enhancer_map, by = 'id') %>%
  filter(gene %in% gene_list$Gene) %>%
  select(-one_of(names(enhancer_map)[c(1, 3, 4, 5, 6)]))

df_atac_enhancer <- atac_enhancer %>%
  pivot_longer(cols = starts_with("TCGA"), names_to = "Sample", values_to = "Value") %>%
  group_by(gene, Sample) %>%
  summarise(Values = list(Value)) %>%
  pivot_wider(names_from = Sample, values_from = Values)

fwrite(df_atac_enhancer, file = "atac_lihc_enhancer.csv")

### Download and process miRNA data
paraCohort <- 'TCGA Kidney Clear Cell Carcinoma'
paraDatasets <- 'miRNA_HiSeq_gene'

mirna_TCGA <- XenaGenerate(subset = XenaHostNames == 'tcgaHub') %>%
  XenaFilter(filterCohorts = paraCohort) %>%
  XenaFilter(filterDatasets = paraDatasets)

XenaQuery(mirna_TCGA) %>% XenaDownload(destdir = 'KIRC')

mirna <- read.csv('KIRC/mirna_kirc', sep = '\t')
selected_samples <- sample(colnames(mirna)[-1], 200)
mirna <- cbind(mirna$sample, mirna[, selected_samples])
names(mirna)[1] <- 'sample'
Accessions <- mirna$sample

# Convert miRNA accessions to names
version <- checkMiRNAVersion(miRNATest$miRNA_Name, verbose = TRUE)
result <- miRNA_AccessionToName(Accessions, targetVersion = "v18")
mirna$sample <- result$TargetName

# Filter relevant miRNA-target pairs
data("mirtarbasegene")
mirna_gene <- subset(mirtarbasegene, Target %in% gene_list$Gene)
mirna <- mirna[, -ncol(mirna)]
names(mirna)[1] <- 'miRNA'
mirna_subset <- subset(mirna, miRNA %in% mirna_gene$miRNA)

mirna_targets <- left_join(mirna_subset, mirna_gene, by = "miRNA")
unique_mirnas <- unique(mirna_targets$miRNA)
all_targets_mirnas <- expand.grid(Target = unique(mirna_gene$Target), miRNA = unique_mirnas)

df_mirna_final <- mirna_targets %>%
  pivot_longer(cols = starts_with("TCGA"), names_to = "Sample", values_to = "Value") %>%
  merge(all_targets_mirnas, by = c("Target", "miRNA"), all = TRUE) %>%
  pivot_wider(names_from = Sample, values_from = Value) %>%
  pivot_longer(cols = starts_with("TCGA"), names_to = "Sample", values_to = "Value") %>%
  mutate(across(starts_with("Value"), ~ replace_na(., 0))) %>%
  group_by(Target, Sample) %>%
  summarise(Values = list(Value)) %>%
  pivot_wider(names_from = Sample, values_from = Values)

fwrite(df_mirna_final, file = "mirna_kirc_df.csv")

### Add logFC and GS annotation
logfc <- read.csv('LIHC/01_DEGsTreat_lihc.csv')
logfc_subset <- subset(logfc, logfc$X %in% gene_list$Gene) %>%
  select(X, logFC) %>% rename(Gene = X)
write.csv(logfc_subset, 'logFC_lihc.csv')

wgcna_gs <- read.csv('LIHC/02_GSandMM_lihc.csv')
wgcna_gs_subset <- subset(wgcna_gs, geneSymbol %in% gene_list$Gene) %>%
  select(geneSymbol, GS_Tumor) %>% rename(Gene = geneSymbol, GS = GS_Tumor)
write.csv(wgcna_gs_subset, 'GS_lihc.csv')
