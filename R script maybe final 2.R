#INSTALL THESE PACKAGE============================================================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("recount3")
BiocManager::install("ComplexHeatmap")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("factoextra")
install.packages("ggpubr")
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

#Creating the sample metadata==================================================================
sample_names <- colnames(indata[, -1]) # Exclude the first column (Gene)
sample.meta.data <- data.frame(
  SRR_ID = sample_names,
  Tissue = sapply(strsplit(sample_names, "_"), `[`, 2), # Extract tissue type (BM or CNS)
  Patient = sapply(strsplit(sample_names, "_"), `[`, 1) # Extract patient ID
)

#Creating the gene metadata====================================================================
gene.meta.data <- data.frame(
  Gene_ID = indata$Gene,
  Gene_Name = indata$Gene) # Using Gene column; replace if another name is available

#Convert data to a count matrix=================================================================
counts_matrix <- as.matrix(indata[, -1]) # Exclude the gene column
rownames(counts_matrix) <- indata$Gene  # Set gene names as rownames

#NORMALIZATION AND FILTERING FOR FILTERED AND NON-FILTERED DATA=================================
counts_matrix <- as.matrix(indata[, -1]) #do not include gene column
rownames(counts_matrix) <- indata$Gene  #set gene names as rownames
dgelist <- DGEList(counts = counts_matrix) #Creates a DGEList object

##Filtering---------------------------------------------------------------------------------
index.keep.expr <- filterByExpr(dgelist, group = sample.meta.data$Tissue)
###Filter dgelist using filtering index, assign it to a new object called dgelist.filtered---
dgelist.filtered <- dgelist[index.keep.expr, , keep.lib.sizes = FALSE]
dim(dgelist.filtered)

###Normalization and filter on filtered list------------------------------------------------
dgelist.filtered.norm <- calcNormFactors(dgelist.filtered, method = "TMM")
dgelist.filtered.norm$samples$norm.factors#view normalization factors

##Calculating LCPM of new filtered list-----------------------------------------------------
lcpm.filtered <- cpm(dgelist.filtered, log = TRUE)

###create a long (tidy) format data frame for LogCPM values--------------------------------------
long.lcpm.filtered <- as.data.frame(lcpm.filtered) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogPCM", names_to = "SRR_ID")

###Create a density plot of the sample distribution-----------------------------------------
ggplot(long.lcpm.filtered) + geom_density(aes(LogPCM, color = SRR_ID)) +
  labs(title = "Log CPM Distribution for Filtered Genes",
       x = "Log CPM",
       y = "Density")

##Combining CPM and LogCPM into a single dataframe================================================
###Calculate CPM and Log2 CPM values from the dgelist.filtered.norm object-----------------
cpm.dgelist.filtered.norm <- cpm(dgelist.filtered.norm)
lcpm.dgelist.filtered.norm <- cpm(dgelist.filtered.norm, log = TRUE)

###Create long format dataframe called df.plotting containing these CPM values-------------
long.cpm.dgelist.filtered.norm <- as.data.frame(cpm.dgelist.filtered.norm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "CPM", names_to = "SRR_ID")

###Create long format dataframe called df.plotting containing these LCPM values------------
long.lcpm.dgelist.filtered.norm <- as.data.frame(lcpm.dgelist.filtered.norm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID")

###Creting a new table that included sample and gene metadata-------------------------------
df.plotting <- full_join(long.cpm.dgelist.filtered.norm, long.lcpm.dgelist.filtered.norm)
df.plotting <- left_join(df.plotting, sample.meta.data, by = c("SRR_ID" = "SRR_ID")) #Joining sample meta data
df.plotting <- left_join(df.plotting, gene.meta.data, by = c("Gene_ID" = "Gene_ID"))

#PCA===============================================================================================================
##Scree plot--------------------------------------------------------------------------------
head(lcpm.dgelist.filtered.norm)
SD <- apply(lcpm.dgelist.filtered.norm, 1 ,sd) #Calculate the SD of the LCPM
threshold <- quantile(SD, 0.9) #threshold for top 10% of SD
top_10_percent <- lcpm.dgelist.filtered.norm[SD >= threshold, ]
##PCA on top 10%----------------------------------------------------------------------------
DS1.svd <- top_10_percent |> 
  t() |> 
  prcomp(scale = FALSE) # PCA using prcomp()
summary(DS1.svd)

fviz_eig(DS1.svd, addlabels = TRUE) + 
  theme_pubr(base_size = 9)

fviz_pca_ind(DS1.svd, label = "none", addEllipses = T, invisible = "quali") +
  theme_pubr(base_size = 9)

##PCA Plot-------------------------------------------------------------------------------------
fviz_pca_ind(DS1.svd, repel = FALSE)



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
top_30_DEGs_Vol <- bind_rows(top_15_BM, top_15_CNS)

###Creating plot--------------------------------------------------------------------------
ggplot() +
  #Plotting significant genes (FDR < 0.05) in blue
  geom_point(
    data = filter(DE, FDR < 0.05),
    aes(x = logFC, y = -log10(PValue)),
    color = "blue", alpha = 0.5, size = 1.5 
  ) +
  #Plotting non-significant genes (FDR > 0.05) in grey
  geom_point(
    data = filter(DE, FDR >= 0.05),
    aes(x = logFC, y = -log10(PValue)),
    color = "grey", alpha= 0.5, size = 1.5
  ) +
  #labelling top 15 DEGs for BM and CNS, colored by Tissue
  geom_text_repel(
    data = top_30_DEGs_Vol,
    aes(x = logFC, y = -log10(PValue), label= Gene_Name, color = Tissue),
    size = 3
  ) +
  #Adding Dashed lines for the thresholds
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  #Customize the plot's appearance
  theme_bw() +
  scale_color_manual(values = c("BM" = "red", "CNS" = "green")) +
  labs(
    title = "Volcano Plot highlighting Top 15 DEGs for BM and CNS",
    subtitle = "Genes colored by tissue type",
    x = "Log2 Fold Change",
    y = "-Log10(PValue)",
    color = "Tissue Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12), 
    legend.position = "top"
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

##Heatmap------------------------------------------------------------------------------------------
###Plot a heatmap of gene cluster 3 values WITHOUT clustering rows or columns, do not plot gene names--
Heatmap(
  z.scaled.genes.cluster3, 
  cluster_rows= FALSE, 
  cluster_columns = FALSE, 
  show_row_names = FALSE
)

###Plot the same data again, clustering columns and rows with the default linkage methods for sample and genes, DO NOT plot gene names---
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

###Subset CPM matrix for top DEGs
cpm_top_DEGs <- cpm.dgelist.filtered.norm[rownames(cpm)]