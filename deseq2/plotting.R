############################################################
# Ferret SARS-CoV-2 RNA-Seq — DESeq2 Analysis Pipeline
# Final_Project2025_2
############################################################

library(tidyverse)
library(DESeq2)
library(pheatmap)

# 1. Load metadata
meta <- read.csv("SraRunTable_fix.csv") %>%
  select(Run, age_group, time, treatment, tissue)
meta <- meta %>% filter(Run != "SRR17047793")

# Check first rows
head(meta)

# 2. Load HTSeq count files
count_dir <- "counts"
sampleFiles <- list.files(count_dir, pattern="*.txt")

# Remove the ".counts.txt" part
sampleNames <- str_remove(sampleFiles, ".counts.txt")

# Check that filenames match metadata Run IDs
stopifnot(all(sampleNames %in% meta$Run))

# Build sampleTable for DESeq2
sampleTable <- data.frame(
  sampleName = meta$Run,
  fileName   = paste0(meta$Run, ".counts.txt"),
  time       = factor(meta$time),
  treatment  = factor(meta$treatment),
  age_group  = factor(meta$age_group)
)

head(sampleTable)

# 3. Construct DESeq2 dataset
#   Design: treatment effect (SARS-CoV-2 vs PBS)

dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory   = count_dir,
  design      = ~ treatment
)

# Filter out genes with very low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

# 4. Results for SARS-CoV-2 vs PBS

res <- results(dds, contrast=c("treatment","SARS-CoV-2","PBS"))
summary(res)

# Order by smallest padj
res_ordered <- res[order(res$padj), ]
head(res_ordered)

# Save full results table
write.csv(as.data.frame(res_ordered),
          "ferret_DESeq2_results.csv",
          row.names=TRUE)

# 5. Variance-stabilized data for PCA & heatmap

vsd <- vst(dds)

# PCA plot
pca_data <- plotPCA(vsd, intgroup=c("treatment","time"), returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color=treatment, shape=time)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

ggplot(pca_data, aes(PC1, PC2, color=treatment, shape=time)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Ferret Lung RNA-Seq Samples") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))


ggsave("PCA_ferret.png", width=6, height=5)

# 6. Heatmap of top 50 most variable genes

# Find most variable genes
rv <- rowVars(assay(vsd))
top50 <- order(rv, decreasing=TRUE)[1:50]

heat_anno <- data.frame(
  treatment = colData(dds)$treatment,
  time      = colData(dds)$time
)
rownames(heat_anno) <- colnames(dds)

pheatmap(
  assay(vsd)[top50,],
  annotation_col = heat_anno,
  show_rownames = FALSE,
  fontsize_row = 6,
  fontsize_col = 8,
  main="Top 50 Most Variable Genes"
)

ggsave("heatmap_top50.png", width=6, height=5)

# 7. Volcano plot

volcano <- as.data.frame(res)
volcano$gene <- rownames(volcano)

ggplot(volcano, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha=0.6) +
  scale_color_manual(values=c("grey","red")) +
  theme_minimal() +
  labs(title="Volcano Plot: SARS-CoV-2 vs PBS")

ggsave("volcano.png", width=6, height=5)

# 8. MA plot

plotMA(res, ylim=c(-6,6))
ggsave("MA_plot.png", width=6, height=5)

############################################################
# END OF SCRIPT

# Gene Enrichment and Pathway Ontology

# Install with bioconductor
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Loading ferret-human ortholog mapping
orth <- read.table("ferret_human_orthologs.tsv.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Keeping relevant columns 
orth_clean <- orth[, c("Gene.stable.ID",
                       "Human.gene.stable.ID",
                       "Human.gene.name")]

colnames(orth_clean) <- c("ferret_id", "human_ensembl", "human_symbol")

head(orth_clean)
summary(orth_clean)

# Merge DE results with orthologs
res_annot <- merge(res_df, orth_clean,
                   by.x = "gene_id", by.y = "ferret_id",
                   all.x = TRUE)

head(res_annot)
dim(res_annot)

# Add a cleaner confidence column
orth$confidence <- as.numeric(orth$Human.orthology.confidence..0.low..1.high.)

# Sort by ferret ID, then by confidence (high first)
orth_sorted <- orth[order(orth$Gene.stable.ID, -orth$confidence), ]

# Deduplicate by ferret gene ID
orth_unique_full <- orth_sorted[!duplicated(orth_sorted$Gene.stable.ID), ]

# Check size
nrow(orth_unique_full)
length(unique(orth_unique_full$Gene.stable.ID))

orth_unique <- orth_unique_full[, c("Gene.stable.ID",
                                    "Human.gene.stable.ID",
                                    "Human.gene.name")]

colnames(orth_unique) <- c("ferret_id", "human_ensembl", "human_symbol")

head(orth_unique)
nrow(orth_unique)

# Merge with unique orthologs
res_annot <- merge(res_df, orth_unique,
                   by.x = "gene_id", by.y = "ferret_id",
                   all.x = TRUE)

dim(res_annot)
sum(is.na(res_annot$human_ensembl))
sum(!is.na(res_annot$human_ensembl))

# Significant DE genes
sig_res <- subset(res_annot, padj < 0.05 & abs(log2FoldChange) > 1)

# DE gene list (human Ensembl IDs)
de_genes <- unique(sig_res$human_ensembl)

# Background gene list (all human orthologs)
bg_genes <- unique(res_annot$human_ensembl)

length(de_genes)
length(bg_genes)

# GO enrichment
ego_bp <- enrichGO(
  gene          = de_genes,
  universe      = bg_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

head(as.data.frame(ego_bp))
dotplot(ego_bp, showCategory = 15, title = "GO Biological Process Enrichment")

# Converting Ensembl IDs to Entrez IDs
gene_entrez <- bitr(de_genes,
                    fromType = "ENSEMBL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

bg_entrez <- bitr(bg_genes,
                  fromType = "ENSEMBL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

de_entrez <- unique(gene_entrez$ENTREZID)
bg_entrez <- unique(bg_entrez$ENTREZID)

# Reactome pathway enrichment
reactome_res <- enrichPathway(
  gene          = de_entrez,
  universe      = bg_entrez,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

head(as.data.frame(reactome_res))

dotplot(reactome_res, showCategory = 15,
        title = "Reactome Pathway Enrichment")


