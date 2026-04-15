rm(list=ls())

library(Matrix)
### HWL_Multilayer network modeling uncovers crucial multi-scale regulatory axes 
### in lung squamous cell carcinoma

### Read example data
setwd(".../DMBN_example")
miRNA_N<-read.csv("Example_data/miRNA_Normal.csv",row.names = 1)
miRNA_D<-read.csv("Example_data/miRNA_Diseased.csv",row.names = 1)
RNA_N<-read.csv("Example_data/RNA_Normal.csv",row.names = 1)
RNA_D<-read.csv("Example_data/RNA_Diseased.csv",row.names = 1)
Protein_N<-read.csv("Example_data/Protein_Normal.csv",row.names = 1)
Protein_D<-read.csv("Example_data/Protein_Diseased.csv",row.names = 1)

### Construct correlation matrix
miRNA_N<-as.matrix(t(miRNA_N))
miRNA_D<-as.matrix(t(miRNA_D))

RNA_N<-as.matrix(t(RNA_N))
RNA_D<-as.matrix(t(RNA_D))

Protein_N<-as.matrix(t(Protein_N))
Protein_D<-as.matrix(t(Protein_D))

miRNA_N_cor <- cor(miRNA_N, use = "all.obs",method="spearman")
miRNA_D_cor <- cor(miRNA_D, use = "all.obs",method="spearman")

RNA_N_cor <- cor(RNA_N, use = "all.obs",method="spearman")
RNA_D_cor <- cor(RNA_D, use = "all.obs",method="spearman")

Protein_N_cor <- cor(Protein_N, use = "all.obs",method="spearman")
Protein_D_cor <- cor(Protein_D, use = "all.obs",method="spearman")

### Construct adjacency matrix
### For ease of demonstration, a threshold of 0.7 is chosen, and threshold 
### discrimination is required before selecting the top r% of correlation coefficients

miRNA_N_adj <- abs(miRNA_N_cor) > 0.7
miRNA_D_adj <- abs(miRNA_D_cor) > 0.7
diff_miRNA_adj<-abs(miRNA_N_adj - miRNA_D_adj)

RNA_N_adj <- abs(RNA_N_cor) > 0.7
RNA_D_adj <- abs(RNA_D_cor) > 0.7
diff_RNA_adj<-abs(RNA_N_adj - RNA_D_adj)

Protein_N_adj <- abs(Protein_N_cor) > 0.7
Protein_D_adj <- abs(Protein_D_cor) > 0.7
diff_Protein_adj<-abs(Protein_N_adj - Protein_D_adj)

### Construct DMBN  ### add intra_edge
DMBN <- as.matrix(bdiag(diff_miRNA_adj, diff_RNA_adj, diff_Protein_adj))

### rename colnames rownames
rownames(DMBN) <- c(rownames(diff_miRNA_adj), rownames(diff_RNA_adj), rownames(diff_Protein_adj))
colnames(DMBN) <- c(colnames(diff_miRNA_adj), colnames(diff_RNA_adj), colnames(diff_Protein_adj))

### Read inter-layer edge data
miRNA_to_RNA<-read.csv("Example_data/miRNA_to_RNA_edge.csv")
RNA_to_Protein<-read.csv("Example_data/RNA_to_Protein_edge.csv")

### add inter_edge
add_edges <- function(mat, edge_df) {
  src <- edge_df$Source
  tgt <- edge_df$Target
  row_idx <- match(src, rownames(mat))
  col_idx <- match(tgt, colnames(mat))
  mat[cbind(row_idx, col_idx)] <- 1
  return(mat)
}

DMBN <- add_edges(DMBN, miRNA_to_RNA)
DMBN <- add_edges(DMBN, RNA_to_Protein)

#write.csv(DMBN, file="Example_data/OutPut/DMBN.csv", row.names=T)



