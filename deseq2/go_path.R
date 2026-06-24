library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

orth <- read.table("ferret_human_orthologs.tsv.txt",
                   +                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

orth_clean <- orth[, c("Gene.stable.ID",
                       "Human.gene.stable.ID",
                       "Human.gene.name")]

colnames(orth_clean) <- c("ferret_id", "human_ensembl", "human_symbol")

head(orth_clean)
summary(orth_clean)

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

res_annot <- merge(res_df, orth_unique,
                   by.x = "gene_id", by.y = "ferret_id",
                   all.x = TRUE)

dim(res_annot)
sum(is.na(res_annot$human_ensembl))
sum(!is.na(res_annot$human_ensembl))

sig_res <- subset(res_annot, padj < 0.05 & abs(log2FoldChange) > 1)

# DE gene list (human Ensembl IDs)
de_genes <- unique(sig_res$human_ensembl)

# Background gene list (all human orthologs)
bg_genes <- unique(res_annot$human_ensembl)

length(de_genes)
length(bg_genes)


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
