# Cell-type and disease enrichment of Cushing genes
setwd("C:/Users/dkeo/surfdrive/cushing_imaging_genetics/cushing_imaging_genetics")
options(stringsAsFactors = FALSE)

library(dplyr)
library(DescTools)
library(ggplot2)
source("functions.R")

# Make output folder
dir.create("output")

########## Load DEGs ########## 

# degs <- read.delim("labelsvulcanoplotCD")
degs <- read.csv("../Final_genes_Bauduin.csv")

# Split into down- and upregulated genes and list entrez IDs
deglist <- list(
  downregulated = as.character(degs[degs$meanofy < 0, "entrez_id"]),
  upregulated = as.character(degs[degs$meanofy > 0, "entrez_id"])
)

total_genes <- 19992

########## Volcano plot ##########

# df <- degs[, c("Gene_symbol", "pvalue", "pvalueadj", "meanofy")]
# df$logp <- -log10(df$pvalueadj)
# ymax <- max(df$logp)
# xmax <- max(abs(df$meanofy))
# df$info <- ifelse(abs(df$meanofy) > 1 & df$pvalueadj < 0.05, '1', '0')
# df$info <- as.factor(df$info)
# df$label <- ""
# # df[mg_genes, "label"] <- entrezId2Name(mg_genes)
# ggplot(df, aes(meanofy, logp, colour = info))+#, label = label)) +
#   geom_point(size = .1) +
#   # geom_text_repel(force = 5, colour = "black", size = 2.5, nudge_y = 0.1, 
#                   # fontface = "italic", segment.size = .1) +
#   scale_colour_manual(values = c("#999999", "#56B4E9")) + # color blind friendly 
#   labs(x = "FC", y = expression('-log'[10]*' '*italic('P')*'-value')) +
#   scale_y_continuous(limits = c(0, ymax)) +
#   scale_x_continuous(limits = c(-xmax, xmax)) +
#   theme_classic() + 
#   theme(legend.position = "none") +
#   ggtitle(r)
# 
# pdf("volcanoplot.pdf", 4, 3)
# vp
# dev.off()

########## Cell-type enrichment ##########

celltype_markers <- readRDS("C:/Users/dkeo/surfdrive/pd_imaging_scn/pd_scn/output/markerlist.rds")
# celltype_markers6 <- celltype_markers[sapply(celltype_markers, length) >= 6]

# Run cell-type enrichment with hyper.test.table
celltype_enrichment <- hypertest.or.table(deglist, celltype_markers, total_genes) # 3d-array: cell-types x measures x deg-type
apply(celltype_enrichment, 3, function(x){
  x[x[, "bh"] < 0.05, ]
})

# Print results
lapply(dimnames(celltype_enrichment)[[3]], function(n){
  x <- celltype_enrichment[,,n]
  x <- cbind('cell-type' = rownames(x), x)
  x <- x[order(x[, "bh"]), ]
  write.table(x, file = paste0("output/celltype_enrichment_", n, ".txt"), quote = FALSE, row.names = FALSE)
})

# # Check any overlap of DEGs and marker genes
# ct_degs <- lapply(celltype_markers, function(l1){
#   lapply(deglist, function(l2){
#     intersect(l1, l2)
#   })
# })
# 
# ct <- sapply(ct_degs, function(m){
#   !any(sapply(m, function(n){
#     length(n)
#   }) == 0)
# })
# ct_degs <- ct_degs[ct]

########## Disease enrichment ##########
disease_markers <- readRDS("C:/Users/dkeo/surfdrive/pd_imaging_scn/pd_scn/output/disease_genes.rds")

disease_enrichment <- hypertest.or.table(deglist, disease_markers, total_genes)
apply(disease_enrichment, 3, function(x){
  x[x[, "bh"] < 0.05, ]
})

# Print results
lapply(dimnames(disease_enrichment)[[3]], function(n){
  x <- disease_enrichment[,,n]
  x <- cbind('disease' = rownames(x), x)
  x <- x[order(x[, "bh"]), ]
  write.table(x, file = paste0("output/disease_enrichment_", n, ".txt"), quote = FALSE, row.names = FALSE)
})

