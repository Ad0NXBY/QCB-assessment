=if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("recount3")
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

#DESCRIPTION OF THE DATASET===================================================================================
indata <- read.delim("B-ALL.txt(1)/B-ALL.txt", sep = "\t")#Read in the data given by Matt for the assessment and assign this to an object called indata.


#Check the number of genes
dimsension <- dim(indata)
nrow(indata)
 
#Extract samples for the dataset that catergorizes the samples into bone marrow and CNS
sample_names <- colnames(indata[-1])
split_sample_names <- str_split_fixed(sample_names, "_", 2)

tissue <- split_sample_names[,2] #Assign this vector to object called tissue


#Extract the part of the sample name that categorizes the patient
patient_number <- str_sub(sample_names,9,9)

#What are the total number of counts per sample?
sample_sums <- colSums(indata[ ,-1])

#How many genes have all zeros for the bone marrow samples or the cerebrospinal fluid samples?
counts_per_gene_sum <- data.frame("Gene" = indata[ ,1],
                                  "BM" = rowSums(indata[ , c(2,4,6)]),
                                  "CNS"= rowSums(indata[ , c(3,5,7)]))

#Try to make use of the pipe to pipe output of one command into the next, filter dataframe to find all zeros from BM and CNS sums
filter(counts_per_gene_sum, BM == 0 | CNS == 0) %>% nrow()

#Using the data in the current format, use base R graphics to plot boxplots of the count distribution per sample
boxplot(indata[, -1])

#Which gene/genes correspond to the very highly expressed outlier?
highly_expressed_genes <- which(indata[,-1] > 3000000, arr.ind = TRUE)
indata[highly_expressed_genes[,1],1]
#Ribosomal RNA removed from data
long.indata <- indata %>%
  as.data.frame() %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Counts") %>%
  mutate("Patients" = str_sub(Sample,9,9),
         "Tissue" = str_split_fixed(Sample, "_",2)[,2])

#NORMALIZATION AND FILTERING==================================================================================
counts_matrix <- as.matrix(indata[, -1]) #do not include gene column
rownames(counts_matrix) <- indata$Gene  #set gene names as rownames
se <- SummarizedExperiment(assays = list(counts = counts_matrix))
counts <- assay(se)
dgelist <- DGEList(counts)

cpm <- cpm(dgelist)
lcpm <- cpm(dgelist, log = TRUE)

#create a long (tidy) format data frame for LogCPM values
long.lcpm <- as.data.frame(lcpm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogPCM", names_to = "SRR_ID")


#Create a density plot of the sample distribution
ggplot(long.lcpm) + geom_density(aes(LogPCM, color = SRR_ID))

#Filtering
index.keep.expr <- filterByExpr(dgelist)

#Filter dgelist using filtering index, assign it to a new object called dgelist.filtered
dgelist.filtered <- dgelist[index.keep.expr, , keep.lib.sizes = FALSE]
dim(dgelist.filtered)

#Plot sample distribution
lcpm.filtered <- cpm(dgelist.filtered, log = TRUE)
long.lcpm.filtered <- as.data.frame(lcpm.filtered) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID")
ggplot(long.lcpm.filtered) + geom_density(aes(LogCPM,color = SRR_ID))

#use calcNormFactors() functions to reate new DGEList object called dgelist.filtered.norm
dgelist.filtered.norm <- calcNormFactors(dgelist.filtered, method = "TMM")
dgelist.filtered.norm$samples$norm.factors

#Calculate CPM and Log2 CPM values from the dgelist.filtered.norm object
cpm.dgelist.filtered.norm <- cpm(dgelist.filtered.norm)
lcpm.dgelist.filtered.norm <- cpm(dgelist.filtered.norm, log = TRUE)

#Create long format dataframe called df.plotting containing these CPM values
long.cpm.dgelist.filtered.norm <- as.data.frame(cpm.dgelist.filtered.norm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "CPM", names_to = "SRR_ID")

#Create long format dataframe called df.plotting containing these LCPM values
long.lcpm.dgelist.filtered.norm <- as.data.frame(lcpm.dgelist.filtered.norm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID")

df.plotting <- full_join(long.cpm.dgelist.filtered.norm, long.lcpm.dgelist.filtered.norm)


w#filters for the gene FCGR1A and plots CPM values split by Celltype and coloured by donor
#ggplot(df.plotting |> filter(Gene_ID == "MIR6723")) + geom_point(aes(x = SRR_ID, y = CPM, color = SRR_ID)) + ggtitle("MIR6723") 
#head(df.plotting)

#PCA ANALYSIS=================================================================================================
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

#IDENTIFICATION OF DEGs BETWEEN GROUPS========================================================================
#Using the DESeq2 Package


#HEATMAP OF TOP 15 DEGS IN EACH GROUP=========================================================================

#BIOLOGICAL INTERPRETATION====================================================================================
