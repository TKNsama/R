# File Name: rna-seq_pipeline_v2.R
# Author: Kenan Tan
# Date: 2024-01-08
# Description: Description of the script

# Load necessary libraries
library(DESeq2)
library(tidyverse)
library(PCAtools)
library(pheatmap)
library(clusterProfiler)
library(KEGGREST)
library(org.At.tair.db)

# Set working directory
setwd("U:/code/r_Pipeline/workspace/")

# setting -------------------------------------------------

# Load raw data
count_matrix <- read.csv("counts_matrix.txt", row.names = 1, header = TRUE, sep = "\t")
coldata <- read.csv("ExpDesignMetaData.csv")
tpm_matrix <- read.csv("tpm.txt", row.names = 1, header = TRUE, sep = "\t")

# Set species and comparison groups
the_species <- "Wheat_v1.1"  # Can be Wheat_v1.1, Wheat_v2.1, Barley_v2 or Barley_v3

unique(coldata$condition)
comparisons <- list(
  c("WT_W2.5_Spike", "bs1_W2.5_Spike"), c("WT_W3_Spike", "bs1_W3_Spike"), c("WT_W3.5_Spike", "bs1_W3.5_Spike"),
  c("WT_W2.5_Node", "bs1_W2.5_Node"), c("WT_W3_Node", "bs1_W3_Node"), c("WT_W3.5_Node", "bs1_W3.5_Node"),
  c("WT_W2.5_Spike", "WT_W3_Spike"), c("WT_W3_Spike", "WT_W3.5_Spike"), c("WT_W2.5_Spike", "WT_W3.5_Spike"),
  c("bs1_W2.5_Spike", "bs1_W3_Spike"), c("bs1_W3_Spike", "bs1_W3.5_Spike"), c("bs1_W2.5_Spike", "bs1_W3.5_Spike"),
  c("WT_W2.5_Node", "WT_W3_Node"), c("WT_W3_Node", "WT_W3.5_Node"), c("WT_W2.5_Node", "WT_W3.5_Node"),
  c("bs1_W2.5_Node", "bs1_W3_Node"), c("bs1_W3_Node", "bs1_W3.5_Node"), c("bs1_W2.5_Node", "bs1_W3.5_Node")
)

# setting end ---------------------------------------------

# Define functions ----------------------------------------

# Differential gene extraction function
diff_gene_list <- function(logfc, comparison, name) {
  res <- results(dds, contrast = c("condition", comparison))
  diff <- subset(res, padj < 0.05 & abs(log2FoldChange) > logfc)
  diff_up <- subset(diff, log2FoldChange > logfc)
  diff_down <- subset(diff, log2FoldChange < -logfc)
  
  # Save results
  write.csv(diff, file = paste0("DEGs/",name, "_all.csv"))
  write.csv(diff_up, file = paste0("DEGs/",name, "_up.csv"))
  write.csv(diff_down, file = paste0("DEGs/",name, "_down.csv"))
}

# Volcano plot generation function
generate_volcano_plot <- function(df, filename) {
  df <- as.data.frame(df)
  x_1 <- df$log2FoldChange
  y_1 <- -log10(df$padj)
  
  # Set colors based on significance
  color <- ifelse(x_1 > 1 & y_1 > -log10(0.05), "red", ifelse(x_1 < -1 & y_1 > -log10(0.05), "blue", "gray"))
  
  # Generate volcano plot
  png(filename = paste0("volcano_plot/",filename, ".jpg"), width = 2400, height = 1800, res = 300)
  plot(x_1, y_1, pch = 20, col = color, main = paste("volcano_plot/",filename, "Volcano Plot"), xlab = "log2FC", ylab = "-log10 Padj")
  dev.off()
}

# GO enrichment analysis function
get_GO_enrich <- function(i) {
  go_enrichment_gene <- read_delim(paste0("U:/code/r_Pipeline/workspace/DEGs/", i), col_names = TRUE)[, 1]
  go_enrichment_gene <- as.factor(go_enrichment_gene$...1)
  go_enrichment_gene_tair <- convert_genes(gene_list = go_enrichment_gene, species = the_species)
  GO_enrich_analysis <- enrichGO(
    gene          = go_enrichment_gene_tair,
    keyType = "TAIR",
    OrgDb         = org.At.tair.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE)
  GO_enrich_analysis_data <- as.data.frame(GO_enrich_analysis)
  
  dir.create(paste("GO/", i, sep = ""))
  
  # Dotplot
  tryCatch({
    pdf(paste("GO/", i,"/GOenrich_dotplot.pdf", sep = ""), width = 10, height = 10)
    print(dotplot(GO_enrich_analysis, showCategory = 10))
    dev.off()
  }, error = function(e) { message("Dotplot error: ", e) })
  dev.off()
  
  # Barplot
  tryCatch({
    pdf(paste("GO/", i,"/GOenrich_barplot.pdf", sep = ""), width = 10, height = 10)
    print(barplot(GO_enrich_analysis, showCategory = 10))
    dev.off()
  }, error = function(e) { message("Barplot error: ", e) })
  dev.off()
  
  write.table(GO_enrich_analysis_data, paste("GO/", i,"/GO_enrichment_results.txt", sep = ""), row.names = FALSE, col.names = TRUE, sep = "\t")
  
  GO_2_classification <- GO_enrich_analysis_data[!is.na(GO_enrich_analysis_data$Description), ]
  GO_2_classification_1_2_3 <- rbind(
    GO_2_classification[GO_2_classification$ONTOLOGY == "BP", ][order(GO_2_classification$p.adjust)[1:10], ],
    GO_2_classification[GO_2_classification$ONTOLOGY == "CC", ][order(GO_2_classification$p.adjust)[1:10], ],
    GO_2_classification[GO_2_classification$ONTOLOGY == "MF", ][order(GO_2_classification$p.adjust)[1:10], ]
  )
  
  # Plotting second classification
  tryCatch({
    pdf(paste("GO/", i,"/GOenrich_second_class.pdf", sep = ""), width = 15, height = 10)
    print(ggplot(GO_2_classification_1_2_3, aes(x = reorder(Description, Count), y = Count, fill = -log10(p.adjust))) +
            geom_bar(stat = "identity") +
            coord_flip() +
            scale_fill_gradient(low = "blue", high = "red") +
            facet_grid(Class ~ ., scale = "free") +
            labs(x = "Term", y = "Counts", fill = "FDR(-log10(P.adjust))") +
            theme(axis.title = element_text(size = 15), axis.text.y = element_text(size = 13)))
    dev.off()
  }, error = function(e) { message("Second classification plot error: ", e) })
  dev.off()
}

# DEGs selection 
All_DEGs_id <- function(folder_path, output_file = NULL) {
  # dir exist?
  if (!dir.exists(folder_path)) {
    stop("The specified folder does not exist.")
  }
  # obtain files path
  file_list <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # csv?
  if (length(file_list) == 0) {
    stop("No CSV files found in the folder.")
  }
  
  columns_list <- lapply(file_list, function(file) {
    tryCatch({
      first_column <- read.csv(file, header = TRUE)[, 1] 
      return(as.character(first_column))
    }, error = function(e) {
      message(paste("Error reading file:", file, "-", e$message))
      return(NULL)
    })
  })
  
  columns_vector <- unlist(columns_list)
  unique_columns <- unique(columns_vector)
  result_df <- data.frame(column = unique_columns, stringsAsFactors = FALSE)
  
  # save result
  if (!is.null(output_file)) {
    write.csv(result_df, output_file, row.names = FALSE)
    message("Result saved to ", output_file)
  }
  
  # get result
  return(result_df)
}

# convert to tair
convert_genes <- function(gene_list, species) {
  
  orthofinder_data <- read.delim("data/Orthogroups.tsv", header = TRUE, sep = "\t")
  
  species_columns <- list(
    Arabidopsis = "Arabidopsis_thaliana_v10.0_representative_gene_model_updated",
    Barley_v2 = "Hordeum_vulgare_Morex_HC_v2.0",
    Barley_v3 = "Hordeum_vulgare_Morex_HC_v3.0",
    Rice = "Oryza_sativa_HC_v7.0",
    Maize = "Zea_mays_B73_NAM_HC_v5.0",
    Wheat_v1.1 = "Triticum_aestivum_CS_HC_v1.1_filtered",
    Wheat_v2.1 = "Triticum_aestivum_CS_HC_v2.1_filtered"
  )
  
  if (!species %in% names(species_columns)) {
    stop("Species not found in the orthofinder data mapping. Please check the species name.")
  }
  
  species_column <- species_columns[[species]]
  
  orthofinder_data <- orthofinder_data %>%
    mutate(
      Arabidopsis_genes = str_split(Arabidopsis_thaliana_v10.0_representative_gene_model_updated, ","),
      species_genes = str_split(!!sym(species_column), ",") # 动态列名
    )
  
  matched_groups <- orthofinder_data %>%
    filter(map_lgl(species_genes, ~ any(. %in% gene_list)))
  
  if(nrow(matched_groups) == 0) {
    message("No matching groups found for the provided gene list.")
  }
  
  unique_arabidopsis_genes <- matched_groups %>%
    pull(Arabidopsis_genes) %>%
    unlist() %>%
    unique()
  
  unique_arabidopsis_genes <- gsub("\\..*", "", unique_arabidopsis_genes)
  
  unique_arabidopsis_genes <- unique(unique_arabidopsis_genes)
  
  if(length(unique_arabidopsis_genes) == 0) {
    message("No matching Arabidopsis genes found.")
  }
  
  return(unique_arabidopsis_genes)
}

# Define functions end ------------------------------------

# Basic process: Data filtering and DESeq2 analysis
count_matrix <- count_matrix[rowSums(count_matrix) > 10, ]
count_matrix <- round(count_matrix)
tpm_matrix <- tpm_matrix[rowSums(tpm_matrix) > 1, ]

# DESeq2 project construction
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

# PCA plot
vst_data <- assay(vst(dds))
rownames(coldata) <- coldata$sample
p <- pca(vst_data, metadata = coldata, removeVar = 0.1)
pca_df <- p$rotated
write.csv(pca_df,"pca_df.csv")
screeplot(p)
biplot(p, x = 'PC1', y = 'PC2')

# Heatmap
pheatmap(tpm_matrix, scale = "row", show_rownames = FALSE)

# Differentially expressed genes and volcano plots
dir.create("DEGs")
for (comp in comparisons) {
  diff_gene_list(logfc = 1, comparison = comp, name = paste0(comp[1], "_vs_", comp[2]))
}

dir.create("volcano_plot")
for (comp in comparisons) {
  comp_name <- paste0(comp[1], "_vs_", comp[2])
  generate_volcano_plot(df = results(dds, contrast = c("condition", comp)), filename = comp_name)
}

# GO enrichment analysis
dir.create("GO")
for (i in dir("U:/code/r_Pipeline/workspace/DEGs/")) {
  get_GO_enrich(i)
}

# Customization --------------------------------------------
# select DEGs for further analysis
DEGs_name <- All_DEGs_id("U:/code/r_Pipeline/workspace/DEGs/",output_file = "Intersection_all_DEGs_id.txt")
DEGs_tpm_matrix <- tpm_matrix[unlist(DEGs_name),]
write.csv(DEGs_tpm_matrix,file = "DEGs_tpm_matrix.csv")
DEGs_tpm_matrix <- read.csv("DEGs_tpm_matrix.csv",row.names = 1)

# Venn plot -----
library(VennDiagram)
library(RColorBrewer)
A <- read.csv("DEGs/WT_W2.5_Spike_vs_bs1_W2.5_Spike_all.csv")[,1]
B <- read.csv("DEGs/WT_W3_Spike_vs_bs1_W3_Spike_all.csv")[,1]
C <- read.csv("DEGs/WT_W3.5_Spike_vs_bs1_W3.5_Spike_all.csv")[,1]
venn.diagram(
  x = list(W2.5=A,W3=B,W3.5=C),
  imagetype="tiff",
  filename = "Spike_all.tiff",
  fill = brewer.pal(5, "Set1")[1:3],
)

# upset plot -----
library(UpSetR)
A <- read.csv("DEGs/resW2.5_Spike_all")[,1]
B <- read.csv("DEGs/resW3_Spike_all")[,1]
C <- read.csv("DEGs/resW3.5_Spike_all")[,1]
E <- read.csv("DEGs/resW2.5_Node_all")[,1]
F1 <- read.csv("DEGs/resW3_Node_all")[,1]
G <- read.csv("DEGs/resW3.5_Node_all")[,1]
max_length <- max(length(A), length(B), length(C), length(E), length(F1), length(G))

A <- c(A, rep(NA, max_length - length(A)))
B <- c(B, rep(NA, max_length - length(B)))
C <- c(C, rep(NA, max_length - length(C)))
E <- c(E, rep(NA, max_length - length(E)))
F1 <- c(F1, rep(NA, max_length - length(F1)))
G <- c(G, rep(NA, max_length - length(G)))

# 创建数据框
df <- data.frame(A, B, C, E, F1, G)
colnames(df) <- c("A","B","C","D","E","F")

upset(fromList(df), sets = c("F","E","D","C","B","A"),
      order.by = 'freq', keep.order = TRUE,
      queries = list(list(query = intersects, 
                          params = list("A","B","C","D","E","F"), 
                          active = T)))

# cluster -----
library(ComplexHeatmap)
library(circlize)
library(cluster)

# K value calculation
sse <- c()
for (k in 1:10) {
  kmeans_model <- kmeans(DEGs_tpm_matrix, centers = k)
  sse <- c(sse, kmeans_model$tot.withinss)
}
plot(1:10, sse, type = "b", pch = 19, frame = FALSE, col = "blue",
     xlab = "Number of clusters (K)", ylab = "SSE",
     main = "Elbow Method for Optimal K")
# Silhouette
asw<-numeric(10) 
for(k in 2:10){ 
  asw[[k]]<-pam(DEGs_tpm_matrix,k)$silinfo$avg.width
}
k.best<-which.max(asw)
cat("通过轮廓系数（silhouette）估计的最佳聚类数：",k.best,"\n") 

# type in the best K here!
K_best <- "3"
normalized_DEGs_tpm_matrix <- t(scale(t(DEGs_tpm_matrix)))
pam_cluster <- pam(normalized_DEGs_tpm_matrix, K_best,metric = "euclidean")
table(pam_cluster$cluster) # see sample size in each cluster
pam_clust<-pam_cluster$cluster
pam_clust<-as.data.frame(pam_clust)
colnames(pam_clust)<-"pam_K10"

# cluster_heatmap
col_fun = colorRamp2(c(-3, 0, 3), c("blue","white","red"))
Heatmap(normalized_DEGs_tpm_matrix,cluster_columns = FALSE,
        col = col_fun,
        #column_split = data.frame(rep(c("aWT_Spike","Bsh1_Spike","cWT_Node","dBsh1_Node"), each = 9)),
        use_raster = F,
        row_split = pam_clust,
        show_row_names = T,border = TRUE,show_heatmap_legend = F,row_names_gp = gpar(fontsize = 8))

# cluster_smooth_line

# in addition -----------------------------------------------

