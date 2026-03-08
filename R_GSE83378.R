#################################################
# ANALISIS DIFFERENTIAL GENE EXPRESSION
# DATASET: GSE83378
# CONDITION: Well-watered vs Drought-stressed
#################################################

############################
# 1. Install & Load Package
############################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery","limma","Biobase"))

install.packages(c("ggplot2","pheatmap","umap","dplyr"))

library(GEOquery)
library(limma)
library(Biobase)
library(ggplot2)
library(pheatmap)
library(umap)
library(dplyr)

############################
# 2. Download Dataset GEO
############################

gset <- getGEO("GSE83378", GSEMatrix = TRUE)

if (length(gset) > 1) idx <- 1 else idx <- 1
gset <- gset[[idx]]

############################
# 3. Expression Matrix
############################

ex <- exprs(gset)

############################
# 4. Ambil Metadata Sampel
############################

pheno <- pData(gset)

# cek metadata
colnames(pheno)

############################
# 5. Tentukan Group
############################

group <- ifelse(
  grepl("drought", pheno$title, ignore.case = TRUE),
  "Drought",
  "Control"
)

group <- factor(group)

table(group)

############################
# 6. Design Matrix
############################

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

############################
# 7. LIMMA ANALYSIS
############################

fit <- lmFit(ex, design)

contrast.matrix <- makeContrasts(
  Drought-Control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

############################
# 8. DEG RESULT
############################

topTableResults <- topTable(
  fit2,
  adjust="fdr",
  number=Inf
)

head(topTableResults)

############################
# 9. BOX PLOT QC
############################

group_colors <- as.numeric(group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Distribusi Ekspresi Gen",
  ylab = "Expression Value"
)

legend(
  "topright",
  legend = levels(group),
  fill = unique(group_colors)
)

############################
# 10. DENSITY PLOT
############################

expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(group, each = nrow(ex))
)

ggplot(expr_long, aes(Expression, color=Group)) +
  geom_density(linewidth=1) +
  theme_minimal()

############################
# 11. UMAP CLUSTERING
############################

umap_input <- t(ex)

umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  Group = group
)

ggplot(umap_df, aes(UMAP1, UMAP2, color=Group)) +
  geom_point(size=3) +
  theme_minimal() +
  ggtitle("UMAP Clustering")

############################
# 12. VOLCANO PLOT
############################

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val
)

volcano_data$status <- "NO"

volcano_data$status[
  volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "UP"

volcano_data$status[
  volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

ggplot(volcano_data,
       aes(logFC, -log10(adj.P.Val), color=status)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("blue","grey","red")) +
  theme_minimal() +
  ggtitle("Volcano Plot")

############################
# 13. HEATMAP TOP 50 DEG
############################

topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults,50)

mat_heatmap <- ex[rownames(top50),]

mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap))==0,
]

gene_variance <- apply(mat_heatmap,1,var)

mat_heatmap <- mat_heatmap[
  gene_variance>0,
]

annotation_col <- data.frame(Group=group)
rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale="row",
  annotation_col=annotation_col,
  show_colnames=FALSE,
  fontsize_row=7,
  main="Top 50 Differentially Expressed Genes"
)

############################
# 14. TOP 20 UP & DOWN GENES
############################

deg <- topTableResults %>%
  filter(adj.P.Val < 0.05)

up_genes <- deg %>%
  arrange(desc(logFC)) %>%
  head(20)

down_genes <- deg %>%
  arrange(logFC) %>%
  head(20)

############################
# 15. SIMPAN HASIL
############################

write.csv(topTableResults,
          "GSE83378_DEG_full.csv")

write.csv(up_genes,
          "GSE83378_top20_upregulated.csv")

write.csv(down_genes,
          "GSE83378_top20_downregulated.csv")

message("Analisis selesai.")
