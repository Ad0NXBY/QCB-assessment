# Load necessary libraries
library(tidyverse)
library(edgeR)
library(ggrepel)

# Read the input data
indata <- read.delim("B-ALL.txt", sep = "\t")

# Extract metadata from column names
sample_names <- colnames(indata[, -1]) # Exclude the first column (Gene)
sample.meta.data <- data.frame(
  SRR_ID = sample_names,
  Tissue = sapply(strsplit(sample_names, "_"), `[`, 2), # Extract tissue type (BM or CNS)
  Patient = sapply(strsplit(sample_names, "_"), `[`, 1) # Extract patient ID
)

# Confirm the structure of sample metadata
head(sample.meta.data)

# Create gene metadata
gene.meta.data <- data.frame(
  Gene_ID = indata$Gene,
  Gene_Name = indata$Gene # Using Gene column; replace if another name is available
)

# Confirm the structure of gene metadata
head(gene.meta.data)

# Convert data to a count matrix
counts_matrix <- as.matrix(indata[, -1]) # Exclude the gene column
rownames(counts_matrix) <- indata$Gene  # Set gene names as rownames

# Create a DGEList object
dgelist <- DGEList(counts = counts_matrix)

# Confirm object structure
class(dgelist)
head(dgelist$counts)
head(dgelist$samples)

# Filter genes with low expression
index.keep.expr <- filterByExpr(dgelist, group = sample.meta.data$Tissue)
head(index.keep.expr)

# Filter the DGEList object
dgelist.filtered <- dgelist[index.keep.expr, , keep.lib.sizes = FALSE]

# Check dimensions of the filtered object
dim(dgelist.filtered)

# Normalize library sizes using TMM
dgelist.filtered.norm <- calcNormFactors(dgelist.filtered, method = "TMM")

# View normalization factors
dgelist.filtered.norm$samples$norm.factors

# Calculate log CPM
lcpm.filtered <- cpm(dgelist.filtered, log = TRUE)

# Convert lcpm to long format
long.lcpm.filtered <- as.data.frame(lcpm.filtered) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID")

# Plot density distribution
ggplot(long.lcpm.filtered) +
  geom_density(aes(LogCPM, colour = SRR_ID)) +
  labs(title = "Log CPM Distribution for Filtered Genes",
       x = "Log CPM",
       y = "Density")

# Combine CPM and logCPM into a single dataframe
cpm.filtered.norm <- cpm(dgelist.filtered.norm)
lcpm.filtered.norm <- cpm(dgelist.filtered.norm, log = TRUE)

long.cpm.filtered.norm <- as.data.frame(cpm.filtered.norm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "CPM", names_to = "SRR_ID")

long.lcpm.filtered.norm <- as.data.frame(lcpm.filtered.norm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID")

df.plotting <- full_join(long.cpm.filtered.norm, long.lcpm.filtered.norm)

# Join sample metadata
df.plotting <- left_join(df.plotting, sample.meta.data, by = c("SRR_ID" = "SRR_ID"))

# Join gene metadata
df.plotting <- left_join(df.plotting, gene.meta.data, by = c("Gene_ID" = "Gene_ID"))

# Confirm structure
head(df.plotting)

# Save necessary objects
#save(dgelist.filtered.norm, sample.meta.data, gene.meta.data,
     #cpm.filtered.norm, lcpm.filtered.norm, file = "Workshop1_output.Rdata")

###ADD PCA HERE###

# Check data types
class(sample.meta.data$Tissue) # Should represent cell type
class(sample.meta.data$Patient) # Should represent donor/patient

# Convert to factors
Patient <- factor(sample.meta.data$Patient)  # Donor equivalent
Tissue <- factor(sample.meta.data$Tissue)  # Celltype equivalent

# Create design matrix for Tissue (cell type equivalent)
design <- model.matrix(~Tissue)

# View the design matrix
print(design)

# Estimate dispersion using the design matrix
dgelist.filtered.norm <- estimateDisp(dgelist.filtered.norm, design = design)

# View dispersion estimates
print(dgelist.filtered.norm$common.dispersion)

# Fit quasi-likelihood model
fit <- glmQLFit(dgelist.filtered.norm, design)

# View the coefficients
head(fit$coefficients)

# Perform differential expression analysis (e.g., CNS vs BM)
BM.V.CNS <- glmQLFTest(fit, coef = 2)  # Coefficient for Tissue (Celltype)

# View results of statistical testing
head(BM.V.CNS$table)

# Extract results using topTags
DE <- topTags(BM.V.CNS, n = Inf)

# Convert to a dataframe
DE <- as.data.frame(DE)

# Annotate with gene names
DE <- left_join(
  rownames_to_column(DE, "gene_id"),  # Add gene_id as a column
  gene.meta.data,  # Join with gene metadata
  by = c("gene_id" = "Gene_ID")
)

# View annotated results
head(DE)

#Volcano plot based on -log10 p value

ggplot(DE, aes(x = logFC, y = -log10(Pvalue)) +
  geom_point(alpha = 0.7) +
  labs(title = "Volcano Plot: Log Fold Change vs -Log10(PValue)",
       x = "Log Fold Change",
       y = "-Log10(PValue)") +
  theme_minimal())

ggplot() +
  geom_point(data = filter(DE, FDR < 0.05), aes(x = logFC, y = -log10(PValue))) +
  geom_point(data = filter(DE, FDR >= 0.05), aes(x = logFC, y = -log10(PValue)), colour = "grey") +
  theme_bw() +
  labs(title = "Volcano Plot with FDR Threshold",
       x = "Log Fold Change",
       y = "-Log10(PValue)")

# Install and load ggrepel if not already installed
# install.packages("ggrepel")
library(ggrepel)

# Select top 15 DEGs for BM (negative logFC) and CNS (positive logFC)
top_15_BM <- DE %>%
  filter(FDR < 0.05, logFC < 0) %>%  # Significant genes for BM
  arrange(logFC) %>%                # Smallest logFC (most negative)
  slice_head(n = 15)                # Top 15 genes

top_15_CNS <- DE %>%
  filter(FDR < 0.05, logFC > 0) %>%  # Significant genes for CNS
  arrange(desc(logFC)) %>%           # Largest logFC (most positive)
  slice_head(n = 15)                 # Top 15 genes

# Combine the top 15 DEGs for BM and CNS
top_30_DEGs <- bind_rows(top_15_BM, top_15_CNS)

# Add a "Tissue" column to indicate the group (BM or CNS)
top_15_BM <- top_15_BM %>%
  mutate(Tissue = "BM")  # Label as BM

top_15_CNS <- top_15_CNS %>%
  mutate(Tissue = "CNS")  # Label as CNS

# Combine the top 15 DEGs for BM and CNS
top_30_DEGs <- bind_rows(top_15_BM, top_15_CNS)

# Volcano plot with colored labels
ggplot() +
  # Plot all significant genes (FDR < 0.05)
  geom_point(data = filter(DE, FDR < 0.05), aes(x = logFC, y = -log10(PValue)), colour = "blue") +
  # Plot non-significant genes (FDR >= 0.05)
  geom_point(data = filter(DE, FDR >= 0.05), aes(x = logFC, y = -log10(PValue)), colour = "grey") +
  # Label top 15 DEGs for each group, colored by Tissue type
  geom_text_repel(
    data = top_30_DEGs,
    aes(x = logFC, y = -log10(PValue), label = Gene_Name, colour = Tissue),
    max.overlaps = Inf
  ) +
  # Customize plot appearance
  theme_bw() +
  scale_color_manual(values = c("BM" = "red", "CNS" = "green")) +  # Set colors for BM and CNS
  labs(
    title = "Volcano Plot with Top 15 DEGs for Each Group (Colored by Tissue Type)",
    x = "Log Fold Change",
    y = "-Log10(PValue)",
    colour = "Tissue Type"
  )



###Heatmap
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")

# Load the Rdata object from Workshop1
load("Workshop1_output.Rdata")

# Rename CPM object for brevity
cpm <- cpm.filtered.norm
head(cpm)

# Create sample names combining Tissue and Patient
sample.meta.data <- sample.meta.data |> 
  mutate(sample_name = paste0(Tissue, "_Patient", Patient))

# Rename columns of CPM
colnames(cpm) <- sample.meta.data$sample_name
head(cpm)

# Filter gene metadata to match CPM genes
filtered.gene.meta.data <- left_join(
  data.frame(gene_id = rownames(cpm)), 
  gene.meta.data, 
  by = c("gene_id" = "Gene_ID")
)

# Rename rows in CPM
rownames(cpm) <- filtered.gene.meta.data$Gene_Name
head(rownames(cpm))

# Perform gene-wise Z-score scaling
z.scaled.genes <- t(cpm) %>% scale() %>% t()
head(z.scaled.genes)

# Calculate Euclidean distances between samples
sample.scaled_distances <- dist(t(z.scaled.genes), method = "euclidean")
print(sample.scaled_distances)

# Calculate Euclidean distances between genes
gene_distances <- dist(z.scaled.genes, method = "euclidean")

# Hierarchical clustering of samples
sample.scaled_hclust <- hclust(sample.scaled_distances, method = "complete")
plot(sample.scaled_hclust, main = "Sample Clustering Dendrogram")

# Hierarchical clustering of genes
gene_hclust <- hclust(gene_distances, method = "average")
plot(gene_hclust, labels = FALSE, main = "Gene Clustering Dendrogram")

# Cut the gene dendrogram into 8 clusters
clusters.genes.k8 <- cutree(gene_hclust, k = 8)

# Table of cluster sizes
table(clusters.genes.k8)

# Extract Z-score data for Cluster 3 genes
z.scaled.genes.cluster3 <- z.scaled.genes[clusters.genes.k8 == 3, ]
dim(z.scaled.genes)
dim(z.scaled.genes.cluster3)

# Save row names (gene names) of Cluster 3 to a CSV
write.csv(rownames(z.scaled.genes.cluster3), file = "cluster3_genenames.csv")

Heatmap(
  matrix = z.scaled.genes.cluster3,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE
)

Heatmap(
  matrix = z.scaled.genes.cluster3,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
)


# Identify top 15 DEGs for BM (logFC < 0) and CNS (logFC > 0)
top_15_BM <- DE %>%
  filter(logFC < 0, FDR < 0.05) %>%  # Significant genes for BM
  arrange(logFC) %>%                # Smallest logFC (most negative)
  slice_head(n = 15)                # Top 15 genes

top_15_CNS <- DE %>%
  filter(logFC > 0, FDR < 0.05) %>%  # Significant genes for CNS
  arrange(desc(logFC)) %>%           # Largest logFC (most positive)
  slice_head(n = 15)                 # Top 15 genes

# Combine top DEGs
top_DEGs <- bind_rows(top_15_BM, top_15_CNS)
top_gene_names <- top_DEGs$Gene_Name

# Subset CPM matrix for top DEGs
cpm_top_DEGs <- cpm.filtered.norm[rownames(cpm.filtered.norm) %in% top_gene_names, ]

# Ensure rows are ordered to match the top DEGs
cpm_top_DEGs <- cpm_top_DEGs[match(top_gene_names, rownames(cpm_top_DEGs)), ]

# Perform Z-score scaling (gene-wise normalization)
z.scaled.cpm_top_DEGs <- t(scale(t(cpm_top_DEGs)))

library(ComplexHeatmap)

# Add sample annotations for Tissue type
sample_annotations <- data.frame(
  Tissue = sample.meta.data$Tissue[match(colnames(z.scaled.cpm_top_DEGs), sample.meta.data$sample_name)]
)

# Convert Tissue to a factor for annotation
sample_annotations$Tissue <- factor(sample_annotations$Tissue, levels = c("BM", "CNS"))

# Define annotation colors
annotation_colors <- list(Tissue = c(BM = "red", CNS = "green"))

# Create a HeatmapAnnotation object
ha <- HeatmapAnnotation(
  df = sample_annotations,
  col = annotation_colors
)

# Draw the heatmap
Heatmap(
  matrix = z.scaled.cpm_top_DEGs,
  top_annotation = ha,          # Add sample annotations
  cluster_rows = TRUE,          # Cluster genes
  cluster_columns = TRUE,       # Cluster samples
  show_row_names = TRUE,        # Show gene names
  show_column_names = TRUE,     # Show sample names
  row_title = "Top 15 DEGs (BM and CNS)",
  column_title = "Samples",
  heatmap_legend_param = list(
    title = "Z-Score",
    legend_direction = "horizontal"
  )
)

# Top 15 DEGs for BM (logFC < 0)
top_15_BM <- DE %>%
  filter(logFC < 0, FDR < 0.05) %>%  # BM-specific genes (negative logFC and significant FDR)
  arrange(logFC) %>%                # Sort by smallest logFC (most negative)
  slice_head(n = 15)                # Select top 15 genes

# Save to CSV
write.csv(top_15_BM, file = "Top_15_DEGs_BM.csv", row.names = FALSE)

# Top 15 DEGs for CNS (logFC > 0)
top_15_CNS <- DE %>%
  filter(logFC > 0, FDR < 0.05) %>%  # CNS-specific genes (positive logFC and significant FDR)
  arrange(desc(logFC)) %>%           # Sort by largest logFC (most positive)
  slice_head(n = 15)                 # Select top 15 genes

# Save to CSV
write.csv(top_15_CNS, file = "Top_15_DEGs_CNS.csv", row.names = FALSE)
