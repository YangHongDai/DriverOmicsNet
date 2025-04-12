# Load necessary libraries
library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(dplyr)

# ---------------------------
# Step 1: Download TCGA Clinical Data
# ---------------------------

# Specify cohort and dataset
tcga_cohort <- 'TCGA Liver Cancer'
tcga_dataset <- 'LIHC_clinicalMatrix'

# Generate query for TCGA clinical data
tcga_clinical_query <- XenaGenerate(subset = XenaHostNames == 'tcgaHub') %>%
  XenaFilter(filterCohorts = tcga_cohort) %>%
  XenaFilter(filterDatasets = tcga_dataset)

# Download clinical matrix for LIHC
XenaQuery(tcga_clinical_query) %>%
  XenaDownload(destdir = 'LIHC')


# ---------------------------
# Step 2: Retrieve GTEx Normal Sample IDs
# ---------------------------

# Load phenotype annotation file
gtex_phenotype <- fread('TcgaTargetGTEX_phenotype.txt')
names(gtex_phenotype) <- gsub("_", "", names(gtex_phenotype))  # Clean column names

# Specify GTEx filtering criteria
gtex_study <- 'GTEX'
gtex_primary_site <- 'Liver'
gtex_tissue_pattern <- '^Liver'

# Subset GTEx for liver samples
gtex_filtered <- subset(gtex_phenotype,
                        study == gtex_study &
                          primarysite == gtex_primary_site &
                          grepl(gtex_tissue_pattern, `primary disease or tissue`))


# ---------------------------
# Step 3: Retrieve TCGA Tumor Sample IDs
# ---------------------------

# Load TCGA clinical matrix
tcga_clinical <- fread(file.path('LIHC', tcga_dataset))
names(tcga_clinical) <- gsub("_", "", names(tcga_clinical))  # Clean column names

# Specify TCGA filtering criteria
tcga_sample_type <- 'Primary Tumor'
tcga_primary_site <- 'Liver'
tcga_histology <- 'Hepatocellular Carcinoma'

# Subset TCGA for primary liver cancer tumor samples
tcga_filtered <- subset(tcga_clinical,
                        sampletype == tcga_sample_type &
                          primarysite == tcga_primary_site &
                          grepl(tcga_histology, histologicaltype))


# ---------------------------
# Step 4: Merge Expression Data from GTEx and TCGA
# ---------------------------

# Select columns (sample IDs) for merging
selected_samples <- c(gtex_filtered$sample, tcga_filtered$sampleID, 'sample')

# Load expression matrix and select relevant samples
expression_data <- fread('TcgaTargetGtex_gene_expected_count', select = selected_samples)

# ---------------------------
# Step 5: Subset to Protein-Coding Genes
# ---------------------------

# Load gene annotation and protein-coding gene list
gene_map <- fread('zz_gencode.v23.annotation.csv', select = c(1, 2))  # Columns: id, gene
protein_coding_genes <- fread('zz_gene.protein.coding.csv')

# Merge expression data with gene map
expression_merged <- merge(gene_map, expression_data, by.x = 'id', by.y = 'sample')

# Keep only protein-coding genes
expression_pc <- subset(expression_merged, gene %in% protein_coding_genes$Gene_Symbol)

# Remove duplicated gene symbols
expression_final <- expression_pc[!(duplicated(expression_pc$gene) | duplicated(expression_pc$gene, fromLast = TRUE)), ]

# Save final expression matrix
write.csv(expression_final, '00_ExpectedCnt_lihc.csv', row.names = FALSE)


# ---------------------------
# Step 6: Prepare Clinical Trait File (Tumor vs. Normal)
# ---------------------------

# Prepare clinical annotation for TCGA (Tumor)
tcga_traits <- tcga_filtered[, c('sampleID', 'primarydisease')]
colnames(tcga_traits) <- c('sampleID', 'tissuetype')
tcga_traits$tissuetype <- 'Tumor'

# Prepare clinical annotation for GTEx (Normal)
gtex_traits <- gtex_filtered[, c('sample', 'sampletype')]
colnames(gtex_traits) <- c('sampleID', 'tissuetype')

# Merge clinical annotation and keep only samples in the expression matrix
clinical_traits <- rbind(tcga_traits, gtex_traits)
clinical_traits_filtered <- subset(clinical_traits, sampleID %in% colnames(expression_final))

# Save final clinical trait matrix
write.csv(clinical_traits_filtered, '00_ClinTraits_lihc.csv', row.names = FALSE)
