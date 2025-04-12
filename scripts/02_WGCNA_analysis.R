# WGCNA Analysis for TCGA Liver Cancer Dataset

library(WGCNA)
options(stringsAsFactors = FALSE)

# Step 1: Load expression data
expr_data <- read.csv("01_voomExpr_lihc.csv")
rownames(expr_data) <- expr_data$X
expr_data <- expr_data[, -1]  # Remove first column (gene names already stored as rownames)

# Step 2: Filter expression data to retain only TCGA samples (optional)
# expr_data <- expr_data[, grepl("TCGA", colnames(expr_data))]

# Step 3: Transpose data (genes as columns, samples as rows)
data_expr <- as.data.frame(t(expr_data))

# Step 4: Quality check
qc_result <- goodSamplesGenes(data_expr, verbose = 3)
stopifnot(qc_result$allOK)

# Step 5: Sample clustering to detect outliers
sample_tree <- hclust(dist(data_expr))
plot(sample_tree, main = "Sample Clustering to Detect Outliers", xlab = "", sub = "")
abline(h = 500, col = "red")

# Step 6: Remove outliers
cut <- cutreeStatic(sample_tree, cutHeight = 500, minSize = 10)
keep_samples <- cut == 1
data_expr <- data_expr[keep_samples, ]

# Step 7: Load clinical trait data
trait_data <- read.csv("00_ClinTraits_lihc.csv")
trait_data$sampleID <- gsub("-", ".", trait_data$sampleID)
trait_data$tissuetype <- ifelse(trait_data$tissuetype == "Tumor", 1, 0)

# Step 8: Match samples between expression and traits
t_match <- match(rownames(data_expr), trait_data$sampleID)
traits <- trait_data[t_match, ]
rownames(traits) <- traits$sampleID
traits <- traits[, "tissuetype", drop = FALSE]

# Step 9: Sample dendrogram with trait coloring
trait_colors <- numbers2colors(traits, signed = FALSE)
plotDendroAndColors(hclust(dist(data_expr)), trait_colors,
                    groupLabels = names(traits),
                    main = "Sample Dendrogram and Trait")

# Step 10: Choose soft-thresholding power
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(data_expr, powerVector = powers, networkType = "signed", verbose = 5)

# Step 11: Plot soft-thresholding results
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit", type = "n",
     main = "Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, col = "red")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red")

# Step 12: Adjacency and TOM matrix
soft_power <- 10
adjacency <- adjacency(data_expr, power = soft_power, type = "signed")
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Step 13: Cluster genes and detect modules
gene_tree <- hclust(as.dist(dissTOM), method = "average")
plot(gene_tree, main = "Gene Clustering Dendrogram")

minModuleSize <- 100
dynamicMods <- cutreeDynamic(dendro = gene_tree, distM = dissTOM,
                             method = "hybrid", deepSplit = 4,
                             pamStage = TRUE, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)

plotDendroAndColors(gene_tree, dynamicColors,
                    "Dynamic Tree Cut", dendroLabels = FALSE,
                    addGuide = TRUE, guideHang = 0.05)

# Step 14: Calculate module eigengenes
MEList <- moduleEigengenes(data_expr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Step 15: Correlate modules with traits
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(data_expr))

# Step 16: Plot module-trait relationships
textMatrix <- paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) <- dim(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-Trait Relationships")

# Step 17: Gene significance and module membership
GS <- as.data.frame(cor(data_expr, traits$tissuetype, use = "p"))
GS_p <- as.data.frame(corPvalueStudent(as.matrix(GS), nrow(data_expr)))
names(GS) <- "GS_Tumor"
names(GS_p) <- "p.GS_Tumor"

mod_names <- substring(names(MEs), 3)
MM <- as.data.frame(cor(data_expr, MEs, use = "p"))
MM_p <- as.data.frame(corPvalueStudent(as.matrix(MM), nrow(data_expr)))
names(MM) <- paste0("MM", mod_names)
names(MM_p) <- paste0("p.MM", mod_names)

# Step 18: Merge and export annotation
annot <- read.csv("zz_gencode.v23.annotation.csv", sep = "\t")
probe_ids <- colnames(t(data_expr))
match_index <- match(probe_ids, annot$gene)

Export <- data.frame(
  geneSymbol = probe_ids,
  geneSymbolCheck = annot$gene[match_index],
  ENSG = annot$id[match_index],
  moduleColor = dynamicColors,
  GS,
  GS_p
)

for (mod in 1:ncol(MM)) {
  Export[[paste0("MM.", mod_names[mod])]] <- MM[, mod]
  Export[[paste0("p.MM.", mod_names[mod])]] <- MM_p[, mod]
}

Export <- Export[order(Export$moduleColor, -abs(Export$GS_Tumor)), ]
write.csv(Export, "02_GSandMM_lihc.csv", row.names = FALSE)
