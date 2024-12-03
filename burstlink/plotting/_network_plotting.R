library(umap)
library(ggplot2)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(readr)
library(ggrepel)
library(factoextra)

setwd("/Users/liying/Desktop/BurstLink-main/docs/tutorials/realistic_data")
gene_umap <- function(counts_file){
  set.seed(123)
  expression_counts_matrix_ <- read_csv(counts_file)
  # UMAP
  expression_counts_matrix <- expression_counts_matrix_[ ,3:ncol(expression_counts_matrix_)]
  expression_counts_matrix <- as.data.frame(lapply(expression_counts_matrix, function(x) as.numeric(as.character(x))))
  umap_config <- umap.defaults
  umap_config$n_neighbors <- 2
  umap_config$min_dist <- 0.5
  umap_result <- umap(expression_counts_matrix, config = umap_config)
  umap_df <- as.data.frame(umap_result$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Gene <- rownames(expression_counts_matrix)
  # Elbow to confirm cluster number
  # fviz_nbclust(umap_df[, c("UMAP1", "UMAP2")], kmeans, method = "wss") + ggtitle("Elbow Method for Determining Optimal Clusters")
  k <- 3
  kmeans_result <- kmeans(umap_df[, c("UMAP1", "UMAP2")], centers = k)
  umap_df$Cluster <- as.factor(kmeans_result$cluster)
  top_genes <- umap_df %>% group_by(Cluster) %>% slice_min(order_by = sqrt(UMAP1^2 + UMAP2^2), n = 5) 
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster, label = Gene)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") + 
    theme_minimal() +
    ggtitle("UMAP Visualization of Genes with Clustering") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_text_repel(max.overlaps = 40, size = 3, box.padding = 0.5, point.padding = 0.5)  
  # print(p)
  return (umap_df)
}

network_plot <- function(gene_interactions_file, burst_info_file, degree_data_file, umap_df){
  gene_interacions <- read_csv(gene_interactions_file)
  burst_info_ <- read_csv(burst_info_file)
  burst_info <- data.frame(Gene = umap_df[ ,3], bs = burst_info_[, 3], bf = log10(burst_info_[, 2]))
  network_data <- data.frame(from = gene_interacions[ ,2], to = gene_interacions[ ,3])
  umap_data <- data.frame(Gene = umap_df[ ,3], UMPA1 = umap_df[ ,1], UMPA2 = umap_df[ ,2])
  umap_data <- umap_data %>% left_join(burst_info, by = "Gene")
  graph <- graph_from_data_frame(d = network_data, vertices = umap_data, directed = TRUE)
  degree_data <- data.frame(indegree = degree(graph, mode = "in"), outdegree = degree(graph, mode = "out"))
  # write.csv(degree_data, file = degree_data_file, row.names = FALSE)
  V(graph)$UMAP1 <- umap_df$UMAP1
  V(graph)$UMAP2 <- umap_df$UMAP2
  V(graph)$Size <- umap_data$X1
  V(graph)$Color <- umap_data$X0
  p <- ggraph(graph, layout = "manual", x = V(graph)$UMAP1, y = V(graph)$UMAP2) +
    geom_edge_fan(alpha = 0.025, color = "gray") +
    geom_node_point(aes(size = V(graph)$Size, color = Color), alpha = 1, shape = 16) +
    scale_size_continuous(range = c(3, 8)) + 
    scale_color_viridis_c() + 
    theme_void()
  print(p)
  # return(list(graph, umap_data))
  return ()
}