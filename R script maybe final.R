library(edgeR)
library(tidyverse)

#Read the input data given from Mat
indata <- read.delim(indata <- read.delim("B-ALL.txt(1)/B-ALL.txt", sep = "\t"))

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
se <- SummarizedExperiment(assays = list(counts = counts_matrix))
counts <- assay(se)
dgelist <- DGEList(counts)

cpm <- cpm(dgelist)
lcpm <- cpm(dgelist, log = TRUE)

###create a long (tidy) format data frame for LogCPM values--------------------------------------
long.lcpm <- as.data.frame(lcpm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogPCM", names_to = "SRR_ID")


###Create a density plot of the sample distribution-----------------------------------------
ggplot(long.lcpm) + geom_density(aes(LogPCM, color = SRR_ID)) #To show before normalization (any significance???)

##Filtering---------------------------------------------------------------------------------
index.keep.expr <- filterByExpr(dgelist)

###Filter dgelist using filtering index, assign it to a new object called dgelist.filtered
dgelist.filtered <- dgelist[index.keep.expr, , keep.lib.sizes = FALSE]
dim(dgelist.filtered)

###Normalization and filter on filtered list------------------------------------------------
dgelist.filtered.norm <- calcNormFactors(dgelist.filtered, method = "TMM")
dgelist.filtered.norm$samples$norm.factors

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


##Calculating LCPM of new filtered list-----------------------------------------------------
lcpm.filtered <- cpm(dgelist.filtered, log = TRUE)
##Plot density distribution by sample
ggplot(long.lcpm.dgelist.filtered.norm) +
  geom_density(aes(LogCPM, colour = SRR_ID)) +
  labs(title = "Log CPM Distribution for Filtered Genes", x = "Log CPM", y = "Density")




#PCA=======================
#Scree plot
head(lcpm.dgelist.filtered.norm)
SD <- apply(lcpm.dgelist.filtered.norm, 1 ,sd) #Calculate the SD of the LCPM
threshold <- quantile(SD, 0.9) #threshold for top 10% of SD
top_10_percent <- lcpm.dgelist.filtered.norm[SD >= threshold, ]
#PCA on top 10%
DS1.svd <- top_10_percent |> 
  t() |> 
  prcomp(scale = FALSE) # PCA using prcomp()
summary(DS1.svd)

fviz_eig(DS1.svd, addlabels = TRUE) + 
  theme_pubr(base_size = 9)

fviz_pca_ind(DS1.svd, label = "none", addEllipses = T, invisible = "quali") +
  theme_pubr(base_size = 9)

#pca plot
fviz_pca_ind(DS1.svd, repel = FALSE)


# Step 9: Create Donor and Celltype factors
Patient <- factor(sample.meta.data$Patient) # Assuming "Patient" represents Donor
Tissue <- factor(sample.meta.data$Tissue) # Assuming "Tissue" represents Celltype

# Step 10: Create a design matrix for Tissue (Celltype)
design <- model.matrix(~Tissue)

# Step 11: Estimate dispersion using the design matrix
dgelist.filtered <- estimateDisp(dgelist.filtered, design = design)

# Display dispersion estimates
print(dgelist.filtered$common.dispersion)

# Step 12: Fit a quasi-likelihood model
fit <- glmQLFit(dgelist.filtered, design)

# Display the first few coefficients
head(fit$coefficients)

# Perform differential expression analysis
BM.vs.CNS <- glmQLFTest(fit, coef = 2) # Coefficient for CelltypeCNS

# View the first few rows of the test results
head(BM.vs.CNS$table)

# Extract results using topTags
DE <- topTags(BM.vs.CNS, n = Inf)

# Convert to a dataframe
DE <- as.data.frame(DE)

# Create gene metadata
gene.meta.data <- data.frame(
  Gene_ID = indata$Gene,
  Gene_Name = indata$Gene # Replace if gene names are available elsewhere
)

# Annotate results with gene metadata
DE <- left_join(
  rownames_to_column(DE, "gene_id"), # Add gene_id as a column
  gene.meta.data, # Join with gene metadata
  by = c("gene_id" = "Gene_ID")
)


library(ggplot2)

# View annotated DE table
head(DE)

# Volcano plot: Log2 Fold Change vs PValue
ggplot(DE, aes(x = logFC, y = PValue)) +
  geom_point(alpha = 0.7) +
  labs(title = "Volcano Plot: Log2 Fold Change vs PValue",
       x = "Log2 Fold Change",
       y = "PValue") +
  theme_minimal()

# Volcano plot: Log2 Fold Change vs -log10(PValue)
ggplot(DE, aes(x = logFC, y = -log10(PValue))) +
  geom_point(alpha = 0.7) +
  labs(title = "Volcano Plot: Log2 Fold Change vs -log10(PValue)",
       x = "Log2 Fold Change",
       y = "-Log10(PValue)") +
  theme_minimal()

