# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import rpy2.robjects as ro


def network_visualization(counts_file, gene_interactions_file, burst_info_file, degree_data_file, network_figure):
    """
    Visualization for gene regulatory interactions network embedded with bursting.
    """
    ro.globalenv['counts_file'] = counts_file
    ro.globalenv['gene_interactions_file'] = gene_interactions_file
    ro.globalenv['burst_info_file'] = burst_info_file
    ro.globalenv['degree_data_file'] = degree_data_file
    ro.globalenv['network_figure'] = network_figure
    network_visualization_R = """
    library(umap)
    library(ggplot2)
    library(ggraph)
    library(igraph)
    library(tidyverse)
    library(RColorBrewer)
    library(readr)
    library(ggrepel)
    library(factoextra)

    set.seed(123)
    expression_counts_matrix_ <- read_csv(counts_file)
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
        scale_size_continuous(name = "Burst size", range = c(3, 8), guide = guide_legend(ticks = FALSE, label = FALSE)) + 
        scale_color_viridis_c(name = "Burst frequency", guide = guide_colorbar(ticks = FALSE, label = FALSE))+
        ggtitle("Gene regulatory interactions network embedded with bursting") +
        theme_void() +
        theme(panel.background = element_rect(fill = "white", color = NA), 
              plot.background = element_rect(fill = "white", color = NA),
              plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
    ggsave(network_figure, plot = p, width = 8, height = 6, dpi = 300)
    """
    ro.r(network_visualization_R)
    return

def TFs_interactiontype_network(readfile_name):
    """
    Visualization for gene regulatory interactions type stacked row chart.
    """
    filtered_inference_result = filter_inference_results(readfile_name, [-3.5, 2.5, -2.5, 5.0])
    feedback_matrix__ = filtered_inference_result[:, np.array([0, 1, 23, 22])]
    unique_genename = sorted(list(set(np.asarray(np.vstack([filtered_inference_result[:, 0], filtered_inference_result[:, 1]])).flatten())))     
    feedback_matrix = np.zeros([len(unique_genename), len(unique_genename)]) - 10
    index = 0
    for genename in unique_genename:
        indices_ = np.where(feedback_matrix__[:, 0] == genename)[0]
        indices__ = np.where(np.isin(unique_genename, feedback_matrix__[indices_, 1]))[0]
        feedback_matrix[indices__, index] = 5 * np.array(feedback_matrix__[indices_, 2].astype(float)).flatten()
        feedback_matrix[index, indices__] = 5 * np.array(feedback_matrix__[indices_, 3].astype(float)).flatten()
        feedback_matrix[index, index] = - 10
        index = index + 1
    feedback_matrix_new = feedback_matrix.copy()
    for index in np.arange(len(unique_genename)):
        feedback_matrix_new[index, :] = sorted(feedback_matrix[index, :], reverse=True)
    feedback_matrix_new = np.hstack([np.arange(len(unique_genename)).reshape([len(unique_genename), 1]), feedback_matrix_new])
    sorted_feedback_matrix =  np.matrix(sorted(feedback_matrix_new.tolist(), key=lambda x: tuple(x[1::])))
    sorted_feedback_matrix_ = np.flipud(sorted_feedback_matrix)
    plt.figure(dpi=300, figsize=(2, 1.5))
    sns.heatmap(sorted_feedback_matrix_[:, 1::], annot=False, cmap=plt.cm.Blues, annot_kws={'size': 100}, linewidths=0.0001, linecolor='w', cbar=False)
    plt.xticks([], [])
    plt.yticks([], [])
    plt.xlabel('Regualtory interaction type', fontsize=3)
    plt.ylabel('TF gene', fontsize=3) 
    plt.title('Regualtory interaction type stacked row chart', fontsize=3)
    plt.show()
    return

def scatterplot_burst(readfile_name, thresholds, group, xlim, ylim):
    """
    Visualization for transcriptional bursting scatter plot.
    """
    filtered_inference_result = filter_inference_results(readfile_name, thresholds)
    bs = np.array(np.vstack([filtered_inference_result[:, 13], filtered_inference_result[:, 14]]).astype(float))
    bf = np.array(np.vstack([filtered_inference_result[:, 11], filtered_inference_result[:, 12]]).astype(float))
    cv2 = np.array(np.vstack([filtered_inference_result[:, 15], filtered_inference_result[:, 16]]).astype(float))
    expression_level = bs * bf

    fig, ax = plt.subplots(dpi=300)
    sc = plt.scatter(np.log10(bf), np.log10(bs), s=cv2*1.5, c=np.log10(expression_level), cmap=plt.cm.viridis)
    cbar = plt.colorbar(sc)
    cbar.ax.tick_params(labelsize=3, width=0.35) 
    cbar.outline.set_linewidth(0.5)
    plt.xlabel('log10(bf)', fontsize=5)
    plt.ylabel('log10(bs)', fontsize=5) 
    if (group == True):
        x = np.linspace(min(np.log10(bf)), max(np.log10(bf)), 100) 
        y = 0.3449 - x  
        ax.plot(x, y, color='red', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=4, width=0.35) 
    ax.tick_params(axis='both', which='minor', labelsize=4, width=0.35) 
    plt.xlim(xlim[0], xlim[1]) 
    plt.ylim(ylim[0], ylim[1])   
    fig.set_size_inches(2, 1.5)
    for spine in ax.spines.values(): 
        spine.set_linewidth(0.5)
    plt.show()
    return

def filter_inference_results(readfile_name, thresholds):
    inference_result_ = pd.read_csv(readfile_name)
    inference_result = np.matrix(inference_result_)[:, 1::]
    indices_infer = np.where(inference_result[:, 17].astype(float) == 0)[0]
    for index in np.arange(4):
        indices_ = np.where(inference_result[:, index+4].astype(float) == 1e-4)[0]
        indices__ = np.where(inference_result[:, index+4].astype(float) == 1e4)[0]
        indices___ = list(set(indices_) | set(indices__))
        indices_infer = list(set(indices_infer) | set(indices___))
    for index in np.arange(2):
        indices_ = np.where(np.log10(inference_result[:, index+11].astype(float)) < thresholds[0])[0]
        indices__ = np.where(np.log10(inference_result[:, index+11].astype(float)) > thresholds[1])[0]
        indices___ = list(set(indices_) | set(indices__))
        indices_infer = list(set(indices_infer) | set(indices___))
    for index in np.arange(2):
        indices_ = np.where(np.log10(inference_result[:, index+13].astype(float)) < thresholds[2])[0]
        indices__ = np.where(np.log10(inference_result[:, index+13].astype(float)) > thresholds[3])[0]
        indices___ = list(set(indices_) | set(indices__))
        indices_infer = list(set(indices_infer) | set(indices___))
    filtered_inference_result = np.delete(inference_result, indices_infer, axis = 0)
    return(filtered_inference_result)