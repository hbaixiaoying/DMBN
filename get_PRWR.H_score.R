
library(igraph)
library(Matrix)

### Read example data
setwd(".../DMBN_example")
DMBN<-read.csv("Example_data/Output/DMBN.csv",row.names = 1)
DMBN <- as.matrix(DMBN)

### directed graph constructed from the adjacency matrix DMBN
g <- graph_from_adjacency_matrix(DMBN, mode = "directed", weighted = NULL, diag = FALSE)

### Compute total degree
out_deg <- degree(g, mode = "out")
in_deg  <- degree(g, mode = "in")
total_deg <- out_deg + in_deg

### Compute betweenness centrality
betw<- betweenness(g, directed = TRUE, normalized = FALSE)

### normalization
normalize_l1 <- function(x) {
  if (sum(x) == 0) {
    return(rep(1/length(x), length(x)))
  } else {
    return(x / sum(x))
  }
}

### Personalized seed score
deg_norm <- normalize_l1(total_deg)
betw_norm <- normalize_l1(betw)
P0_unnorm <- deg_norm + betw_norm
P0 <- normalize_l1(P0_unnorm)

### node id
miRNA_nodes <- rownames(DMBN)[1:20]
RNA_nodes   <- rownames(DMBN)[21:65]
Protein_nodes <- rownames(DMBN)[66:100]

### miRNA layer
DMBN_mm <- DMBN[miRNA_nodes, miRNA_nodes]

### RNA layer
DMBN_rr <- DMBN[RNA_nodes, RNA_nodes]

### Protein layer
DMBN_pp <- DMBN[Protein_nodes, Protein_nodes]

### miRNA -> RNA
DMBN_mr <- DMBN[miRNA_nodes, RNA_nodes]

### miRNA -> Protein
DMBN_mp <- DMBN[miRNA_nodes, Protein_nodes]

### RNA -> miRNA
DMBN_rm <- DMBN[RNA_nodes, miRNA_nodes]

### RNA -> Protein
DMBN_rp <- DMBN[RNA_nodes, Protein_nodes]

### Protein -> miRNA
DMBN_pm <- DMBN[Protein_nodes, miRNA_nodes]

### Protein -> RNA
DMBN_pr <- DMBN[Protein_nodes, RNA_nodes]

### Store all blocks in a list
network_blocks <- list(
  mm = DMBN_mm,
  rr = DMBN_rr,
  pp = DMBN_pp,
  mr = DMBN_mr,
  mp = DMBN_mp,
  rm = DMBN_rm,
  rp = DMBN_rp,
  pm = DMBN_pm,
  pr = DMBN_pr
)

### PRWR-H algorithm
multilayer_RWR_strict <- function(adj_blocks, lambda = 0.7) {
  
  # ---------- some information ----------
  n_mirna    <- nrow(adj_blocks$mm)
  n_rna      <- nrow(adj_blocks$rr)
  n_protein  <- nrow(adj_blocks$pp)
  total_nodes <- n_mirna + n_rna + n_protein
  
  all_nodes <- c(
    rownames(adj_blocks$mm),
    rownames(adj_blocks$rr),
    rownames(adj_blocks$pp)
  )
  
  H <- Matrix(0, nrow = total_nodes, ncol = total_nodes, sparse = TRUE)
  rownames(H) <- all_nodes
  colnames(H) <- all_nodes
  
  miRNA_idx   <- 1:n_mirna
  RNA_idx     <- (n_mirna + 1):(n_mirna + n_rna)
  Protein_idx <- (n_mirna + n_rna + 1):total_nodes
  
  # ---------- helper function ----------
  norm_vec <- function(x) {
    s <- sum(x)
    if (s == 0) return(x)
    x / s
  }
  
  # ---------- miRNA layer ----------
  for (i in seq_len(n_mirna)) {
    
    k <- 0
    if (sum(adj_blocks$mr[i, ]) > 0) k <- k + 1
    if (sum(adj_blocks$mp[i, ]) > 0) k <- k + 1
    
    # intra layer
    if (k > 0) {
      H[miRNA_idx[i], miRNA_idx] <-
        (1 - lambda) * norm_vec(adj_blocks$mm[i, ])
    } else {
      H[miRNA_idx[i], miRNA_idx] <-
        norm_vec(adj_blocks$mm[i, ])
    }
    
    # inter layer
    if (k > 0) {
      w <- lambda / k
      if (sum(adj_blocks$mr[i, ]) > 0) {
        H[miRNA_idx[i], RNA_idx] <-
          w * norm_vec(adj_blocks$mr[i, ])
      }
      if (sum(adj_blocks$mp[i, ]) > 0) {
        H[miRNA_idx[i], Protein_idx] <-
          w * norm_vec(adj_blocks$mp[i, ])
      }
    }
  }
  
  # ---------- RNA layer ----------
  for (i in seq_len(n_rna)) {
    
    k <- 0
    if (sum(adj_blocks$rm[i, ]) > 0) k <- k + 1
    if (sum(adj_blocks$rp[i, ]) > 0) k <- k + 1
    
    # intra layer
    if (k > 0) {
      H[RNA_idx[i], RNA_idx] <-
        (1 - lambda) * norm_vec(adj_blocks$rr[i, ])
    } else {
      H[RNA_idx[i], RNA_idx] <-
        norm_vec(adj_blocks$rr[i, ])
    }
    
    # inter layer
    if (k > 0) {
      w <- lambda / k
      if (sum(adj_blocks$rm[i, ]) > 0) {
        H[RNA_idx[i], miRNA_idx] <-
          w * norm_vec(adj_blocks$rm[i, ])
      }
      if (sum(adj_blocks$rp[i, ]) > 0) {
        H[RNA_idx[i], Protein_idx] <-
          w * norm_vec(adj_blocks$rp[i, ])
      }
    }
  }
  
  # ---------- Protein layer ----------
  for (i in seq_len(n_protein)) {
    
    k <- 0
    if (sum(adj_blocks$pm[i, ]) > 0) k <- k + 1
    if (sum(adj_blocks$pr[i, ]) > 0) k <- k + 1
    
    # intra layer
    if (k > 0) {
      H[Protein_idx[i], Protein_idx] <-
        (1 - lambda) * norm_vec(adj_blocks$pp[i, ])
    } else {
      H[Protein_idx[i], Protein_idx] <-
        norm_vec(adj_blocks$pp[i, ])
    }
    
    # inter layer
    if (k > 0) {
      w <- lambda / k
      if (sum(adj_blocks$pm[i, ]) > 0) {
        H[Protein_idx[i], miRNA_idx] <-
          w * norm_vec(adj_blocks$pm[i, ])
      }
      if (sum(adj_blocks$pr[i, ]) > 0) {
        H[Protein_idx[i], RNA_idx] <-
          w * norm_vec(adj_blocks$pr[i, ])
      }
    }
  }
  
  return(H)
}

### iteration update
rwr_iteration <- function(transition_matrix, restart_prob = 0.5, max_iter = 3000, tol = 1e-6) {
  all_nodes <- rownames(transition_matrix)
  n <- length(all_nodes)
  #Personalized seed
  p <- P0
  
  for (i in 1:max_iter) {     
    p_next <- (1 - restart_prob) * (t(transition_matrix) %*% p) + restart_prob * p
    
    if (sum(abs(p_next - p)) < tol) {
      cat("Converged at iteration", i, "\n")
      return(as.vector(p_next))
    }
    p <- p_next
  }
  
  cat("Maximum number of iterations reached:", max_iter, "\n")
  return(as.vector(p))
}

### transition_matrix
transition_matrix <- multilayer_RWR_strict(network_blocks, lambda = 0.7)

### operation RWR-H
PRWRH_scores <- rwr_iteration(transition_matrix,  restart_prob = 0.5)
PRWRH <- data.frame( Node = rownames(transition_matrix), PRWRH_Score = PRWRH_scores )

#write.csv(PRWRH, file="Example_data/OutPut/PRWRH_score.csv", row.names=F)

