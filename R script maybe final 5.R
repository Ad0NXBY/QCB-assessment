#INSTALL THESE PACKAGE============================================================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "recount3", "ComplexHeatmap"))
install.packages(c("tidyverse", "ggplot2", "factoextra", "ggpubr"))

library(edgeR)
library(recount3)
library(tidyverse)
library(ggplot2)
library(factoextra)
library(ggpubr)
library(ComplexHeatmap)
library(ggrepel)

#Read the input data==========================================================================
indata <- read.delim("B-ALL.txt(1)/B-ALL.txt", sep = "\t")

#Creating the sample metadata ================================================================
sample_names <- colnames(indata[, -1]) # Exclude the first column (Gene)
sample.meta.data <- data.frame(
  SRR_ID = sample_names,
  Tissue = sapply(strsplit(sample_names, "_"), `[`, 2), # Extract tissue type (BM or CNS)
  Patient = sapply(strsplit(sample_names, "_"), `[`, 1) # Extract patient ID
)

#Creating the gene metadata===================================================================
gene.meta.data <- data.frame(
  Gene_ID = indata$Gene,
  Gene_Name = indata$Gene # Using Gene column; replace if another name is available
)

##Convert data to a count matrix --------------------------------------------------------------
counts_matrix <- as.matrix(indata[, -1]) # Exclude the gene column
rownames(counts_matrix) <- indata$Gene  # Set gene names as rownames

#Normalization and filtering for filtered and non-filtered data ==============================
dgelist <- DGEList(counts = counts_matrix) # Creates a DGEList object

##Filtering-----------------------------------------------------------------------------------
index.keep.expr <- filterByExpr(dgelist, group = sample.meta.data$Tissue) # FIXED: Correct grouping

##Filter lowly expressed genes ----------------------------------------------------------------
dgelist.filtered <- dgelist[index.keep.expr, , keep.lib.sizes = FALSE]

##Normalization and filtering on filtered list -----------------------------------------------
dgelist.filtered.norm <- calcNormFactors(dgelist.filtered, method = "TMM")

##Calculate LCPM of new filtered list --------------------------------------------------------
lcpm.filtered <- cpm(dgelist.filtered, log = TRUE)

##Create a long (tidy) format data frame for LogCPM values ------------------------------------
long.lcpm.filtered <- as.data.frame(lcpm.filtered) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogPCM", names_to = "SRR_ID")

##Create a density plot of the sample distribution -------------------------------------------
ggplot(long.lcpm.filtered) + geom_density(aes(LogPCM, color = SRR_ID)) +
  labs(title = "Log CPM Distribution for Filtered Genes",
       x = "Log CPM",
       y = "Density")

#Combining CPM and LogCPM into a single dataframe ===========================================
##Calculate CPM and Log2 CPM values from the dgelist.filtered.norm object --------------------
cpm.dgelist.filtered.norm <- cpm(dgelist.filtered.norm)
lcpm.dgelist.filtered.norm <- cpm(dgelist.filtered.norm, log = TRUE)

#Boxplots to show before and after normalization ============================================
##Calculate LogCPM for raw counts (before normalization, use filtered data) ------------------
raw.lcpm <- cpm(dgelist.filtered, log = TRUE) # FIXED: Use filtered data for consistency

##Create a long format dataframe for raw LogCPM values (before normalization)
long.raw.lcpm <- as.data.frame(raw.lcpm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID") %>%
  mutate(Status = "Before Normalization")

##Create a long format dataframe for normalized LogCPM values (after normalization)
long.norm.lcpm <- as.data.frame(lcpm.dgelist.filtered.norm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID") %>%
  mutate(Status = "After Normalization")

##Combine both data frames for plotting -----------------------------------------------------
combined.lcpm <- bind_rows(long.raw.lcpm, long.norm.lcpm)

##Plot the boxplots --------------------------------------------------------------------------
ggplot(combined.lcpm, aes(x = SRR_ID, y = LogCPM, fill = Status)) +
  geom_boxplot(outlier.shape = 0.5, alpha = 0.5) + # Boxplot without outliers
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x-axis labels
  facet_wrap(~ Status, nrow = 2) + # Separate before and after normalization
  labs(
    title = "Comparison of LogCPM Before and After Normalization",
    x = "Sample",
    y = "LogCPM"
  ) +
  theme_minimal()



#PCA===============================================================================================================
##Scree plot--------------------------------------------------------------------------------
DS1.svd <- lcpm.dgelist.filtered.norm |> 
  t() |> 
  prcomp(scale = FALSE) # PCA using prcomp()
summary(DS1.svd)

fviz_eig(DS1.svd, addlabels = TRUE) + 
  theme_pubr(base_size = 9)

##PCA Plot-------------------------------------------------------------------------------------
fviz_pca_ind(DS1.svd, repel = FALSE) + 
  expand_limits(x = 6) +
  labs(title = "B_ALL PCA",
       x = "PC1",
       y = "PC2")


#IDENTIFICATION OF DEGs BETWEEN GROUPS==============================================================================
##Check the data types------------------------------------------------------------------------
class(sample.meta.data$Tissue)
class(sample.meta.data$Patient)

##Create Donor and Cell type factors----------------------------------------------------------
Patient <- factor(sample.meta.data$Patient) # Assuming "Patient" represents Donor
Tissue <- factor(sample.meta.data$Tissue) # Assuming "Tissue" represents Celltype

##Create a design matrix for Tissue (Celltype)-------------------------------------------------
design <- model.matrix(~Tissue)

##Estimate dispersion using the design matrix--------------------------------------------------
dgelist.filtered <- estimateDisp(dgelist.filtered, design = design)

###Display dispersion estimates-------------------------------------------------------------
print(dgelist.filtered$common.dispersion)


##Fit a quasi-likelihood model--------------------------------------------------------------
qlfit <- glmQLFit(dgelist.filtered, design)

###Display the first few coefficients-------------------------------------------------------
head(qlfit$coefficients)

##Perform differential expression analysis-------------------------------------------------
BM.vs.CNS <- glmQLFTest(qlfit, coef = 2) # Coefficient for CelltypeCNS

###View the first few rows of the test results---------------------------------------------
head(BM.vs.CNS$table)

##Extract results using topTags------------------------------------------------------------
DE <- topTags(BM.vs.CNS, n = Inf)

###Convert to a dataframe------------------------------------------------------------------
DE <- as.data.frame(DE)

##Annotate results with gene names from gene.meta.data using left join---------------------
DE <- left_join(
  rownames_to_column(DE, "gene_id"), # Add gene_id as a column
  gene.meta.data, # Join with gene metadata
  by = c("gene_id" = "Gene_ID")
)

###View annotated DE table-----------------------------------------------------------------
head(DE)


##Volcano plots---------------------------------------------------------------------------
###Selecting the top 15 DEgs for BM and CNS respectively----------------------------------
top_15_BM_Vol <- DE %>%
  filter(FDR < 0.05, logFC < 0) %>%  # Significant genes with negative logFC
  arrange(logFC) %>%                 # Smallest logFC (most negative)
  slice_head(n = 15) %>%             # Top 15 genes
  mutate(Tissue = "BM")              # Label as BM

top_15_CNS_Vol <- DE %>%
  filter(FDR < 0.05, logFC > 0) %>%  # Significant genes with positive logFC
  arrange(desc(logFC)) %>%           # Largest logFC (most positive)
  slice_head(n = 15) %>%             # Top 15 genes
  mutate(Tissue = "CNS")             # Label as CNS

###Combining the top 15 DEGs for BM and CNS-----------------------------------------------
top_30_DEGs_Vol <- bind_rows(top_15_BM_Vol, top_15_CNS_Vol)

### Creating plot--------------------------------------------------------------------------
### Creating plot--------------------------------------------------------------------------
ggplot() +
  # Plot non-significant genes (FDR >= 0.05 or -1 < logFC < 1) in grey
  geom_point(
    data = filter(DE, FDR >= 0.05 | (logFC > -1 & logFC < 1)),
    aes(x = logFC, y = -log10(PValue), color = "Non-Significant"),
    alpha = 0.5, size = 1.5
  ) +
  # Plot downregulated significant genes (FDR < 0.05 and logFC <= -1) in orange
  geom_point(
    data = filter(DE, FDR < 0.05 & logFC <= -1),
    aes(x = logFC, y = -log10(PValue), color = "Downregulated"),
    alpha = 0.5, size = 1.5
  ) +
  # Plot upregulated significant genes (FDR < 0.05 and logFC >= 1) in blue
  geom_point(
    data = filter(DE, FDR < 0.05 & logFC >= 1),
    aes(x = logFC, y = -log10(PValue), color = "Upregulated"),
    alpha = 0.5, size = 1.5
  ) +
  # Label top 15 DEGs for BM and CNS, colored by Tissue
  geom_text_repel(
    data = top_30_DEGs_Vol,
    aes(x = logFC, y = -log10(PValue), label = Gene_Name, color = Tissue),
    size = 3, max.overlaps = 20  # Increased max.overlaps to handle overlaps better
  ) +
  # Adding dashed lines for the thresholds
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.02), linetype = "dashed", color = "black", linewidth = 0.8) +
  # Customize the plot's appearance
  theme_bw() +
  scale_color_manual(
    values = c(
      "Upregulated" = "blue", 
      "Downregulated" = "orange", 
      "Non-Significant" = "grey",
      "BM" = "red", 
      "CNS" = "green"
    )
  ) +
  labs(
    title = "Volcano Plot highlighting Top 15 DEGs for BM and CNS",
    subtitle = "Genes colored by tissue type and regulation",
    x = "Log2 Fold Change",
    y = "-Log10(PValue)",
    color = "Gene Category"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",  # Move the legends to the right
    legend.box = "vertical",    # Stack legends vertically
    legend.box.spacing = unit(0.5, "cm")  # Add spacing between the legends
  ) +
  guides(
    color = guide_legend(title = "Gene Category", override.aes = list(size = 3)),  # Gene Category legend
    fill = guide_legend(title = "Tissue Type")  # Tissue Type legend
  )



#HEATMAP OF TOP 15 DEGS IN EACH GROUP==============================================================
##Renamine CPM.dgelist.filtered.norm to cpm for ease-----------------------------------------------
cpm <- cpm.dgelist.filtered.norm

##Create sample names combining Tissue and Patient-------------------------------------------------
sample.meta.data <- sample.meta.data %>%
  mutate(sample_name = paste0(Tissue, "_Patient", Patient))

##Rename columns of CPM----------------------------------------------------------------------------
colnames(cpm) <- sample.meta.data$sample_name
head(cpm)#to check that everything is going correct

###Filter gene metadata to match CPM genes------------------------------------------------------
filtered.gene.meta.data <- left_join(
  data.frame(gene_id = rownames(cpm)),
  gene.meta.data,
  by = c("gene_id" = "Gene_ID")
)

##Rename rows of CPM------------------------------------------------------------------------------
rownames(cpm) <- filtered.gene.meta.data$Gene_Name

##Use Z-Score scaling on the cpm data and assign this to an object called z.scaled.genes----------
z.scaled.genes <- t(cpm) %>% scale () %>% t() 
head(z.scaled.genes)

##Find Euclidean distance for between samples and Hierarchical clustering of samples-------------
sample.scaled_distance <- dist(t(z.scaled.genes), method = "euclidean") 
sample.scale_hclust <- hclust(sample.scaled_distance, method = "complete")
plot(sample.scale_hclust, main = "Sample Clustering Dendrogram")

##Find Euclidean distance for between genes and Hierarchical clustering of genes-----------------
gene_distance <- dist(z.scaled.genes, method = "euclidean")
gene_hclust <- hclust(gene_distance, method = "average")
plot(gene_hclust, labels = FALSE, main = "Gene Clustering Dendrogram")

##Cut the gene dendrogram into 8 clusters--------------------------------------------------------
clusters.gene.k8 <- cutree(gene_hclust, k = 8)
table(clusters.gene.k8) #table of cluster sizes

##Subset the z-score scaled data for the genes that are in the gene cluster 3--------------------
z.scaled.genes.cluster3 <- z.scaled.genes[clusters.gene.k8 ==3, ]

###Check the dimensions of the subsetted data----------------------------------------------------
dim(z.scaled.genes)
dim(z.scaled.genes.cluster3)

##Save row names to a csv file-------------------------------------------------------------------
write.csv(rownames(z.scaled.genes.cluster3), file = "cluster3_genenames.csv")

#HEATMAP==========================================================================================
##Plot a heatmap of gene cluster 3 values WITHOUT clustering rows or columns, do not plot gene names--
Heatmap(
  z.scaled.genes.cluster3, 
  cluster_rows= FALSE, 
  cluster_columns = FALSE, 
  show_row_names = FALSE
)

##Plot the same data again, clustering columns and rows with the default linkage methods for sample and genes, DO NOT plot gene names---
Heatmap(
  z.scaled.genes.cluster3,
  cluster_rows= TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE
)


##Identify the top 15 DEGs for BM(logFC < 0) and CND (logFC > 0) and combine the top DEGs------------
top_15_BM_Heat <- DE %>%
  filter(FDR < 0.05, logFC < 0) %>%  
  arrange(logFC) %>%                
  slice_head(n = 15) 

top_15_CNS_Heat <- DE %>%
  filter(FDR < 0.05, logFC > 0) %>%  
  arrange(desc(logFC)) %>%          
  slice_head(n = 15)

top_DEGs_Heat <- bind_rows(top_15_BM_Heat, top_15_CNS_Heat)
top_gene_names_Heat <- top_DEGs_Heat$Gene_Name

##Subset CPM matrix for top DEGs----------------------------------------------------------------------
lcpm_top_DEGs_Heat <- lcpm.dgelist.filtered.norm[rownames(lcpm.dgelist.filtered.norm) %in% top_gene_names_Heat, ]

##Ensuring rows are ordered to match the top DEGs----------------------------------------------------
lcpm_top_DEGs_Heat <- lcpm_top_DEGs_Heat[match(top_gene_names_Heat, rownames(lcpm_top_DEGs_Heat)), ]

##Perform Z- score scaling (gene-wise normalization)--------------------------------------------------
z.scaled.lcpm_top_DEGs_Heat <- t(scale(t(lcpm_top_DEGs_Heat)))

##Heatmap
Heatmap(
  matrix = z.scaled.lcpm_top_DEGs_Heat,
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


