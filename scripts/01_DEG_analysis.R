# Load required libraries
library(dplyr)
library(limma)
library(edgeR)

# -------------------------------
# Step 1: Back-transform log2(count + 1) values to raw count
# -------------------------------

# Load preprocessed expression matrix
expr_df <- read.csv('00_ExpectedCnt_lihc.csv')

# Remove metadata columns and retain count matrix
expr_counts <- expr_df[, -c(1:3)]
rownames(expr_counts) <- expr_df$gene

# Back-transform to approximate raw count values
expr_counts_bt <- round((2^expr_counts - 1), 0)
write.csv(expr_counts_bt, '01_ExpectedCntBT_lihc.csv')

# -------------------------------
# Step 2: Create DGEList object
# -------------------------------

dge <- DGEList(counts = expr_counts_bt)

# Group samples by origin (TCGA tumor vs. GTEx normal)
sample_ids <- colnames(dge)
group <- substr(sample_ids, 1, 4)
dge$samples$group <- group

# -------------------------------
# Step 3: Quality control and filtering
# -------------------------------

# Log-CPM transformation for diagnostic plots
lcpm <- cpm(dge, log = TRUE)

# Examine library size distribution
cat("Mean Library Size (Million):", mean(dge$samples$lib.size) * 1e-6, "\n")
cat("Median Library Size (Million):", median(dge$samples$lib.size) * 1e-6, "\n")

# Check group sizes
table(dge$samples$group)

# Filter lowly-expressed genes
keep_genes <- filterByExpr(dge, group = group)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

cat("Remaining genes after filtering:", dim(dge)[1], "\n")

# Plot density of filtered logCPM
lcpm <- cpm(dge, log = TRUE)
plot(density(lcpm[, 1]), lwd = 2, ylim = c(0, 0.26), las = 2,
     main = "Filtered Log-CPM Density", xlab = "Log-CPM")

# -------------------------------
# Step 4: Normalize and model design
# -------------------------------

# Normalize library size (TMM method)
dge <- calcNormFactors(dge)

# Create design matrix for voom-limma
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))

# Define contrast: Tumor (TCGA) vs. Normal (GTEX)
contrast_matrix <- makeContrasts(TCGAvsGTEX = TCGA - GTEX, levels = colnames(design))

# -------------------------------
# Step 5: Apply voom transformation
# -------------------------------

v <- voom(dge, design, plot = TRUE)

# Fit linear model and apply contrasts
fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrasts = contrast_matrix)

# Apply empirical Bayes moderation
fit_ebayes <- eBayes(fit)

# Plot mean-variance trend
plotSA(fit_ebayes, main = "Mean-Variance Trend (eBayes)")

# Examine initial DEG result (uncorrected)
summary(decideTests(fit_ebayes))

# -------------------------------
# Step 6: Refine DEG detection using treat (logFC threshold = 1)
# -------------------------------

fit_treat <- treat(fit, lfc = 1)
summary(decideTests(fit_treat))

# Extract DEG table
deg_result <- topTreat(fit_treat, n = Inf)

# Save results
write.csv(deg_result, "01_DEGsTreat_lihc.csv")

# -------------------------------
# Step 7: Export voom-normalized expression matrix for WGCNA
# -------------------------------

voom_expr_matrix <- v$E
write.csv(voom_expr_matrix, '01_voomExpr_lihc.csv')
