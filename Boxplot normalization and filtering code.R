# Read the input data
indata <- read.delim("B-ALL.txt(1)/B-ALL.txt", sep = "\t")

# Creating the sample metadata ================================================================
sample_names <- colnames(indata[, -1]) # Exclude the first column (Gene)
sample.meta.data <- data.frame(
  SRR_ID = sample_names,
  Tissue = sapply(strsplit(sample_names, "_"), `[`, 2), # Extract tissue type (BM or CNS)
  Patient = sapply(strsplit(sample_names, "_"), `[`, 1) # Extract patient ID
)

# Validate metadata to ensure it matches expectations
print(head(sample.meta.data)) # FIXED: Added for validation

# Creating the gene metadata ------------------------------------------------------------------
gene.meta.data <- data.frame(
  Gene_ID = indata$Gene,
  Gene_Name = indata$Gene # Using Gene column; replace if another name is available
)

# Convert data to a count matrix --------------------------------------------------------------
counts_matrix <- as.matrix(indata[, -1]) # Exclude the gene column
rownames(counts_matrix) <- indata$Gene  # Set gene names as rownames

# Normalization and filtering for filtered and non-filtered data ==============================
dgelist <- DGEList(counts = counts_matrix) # Creates a DGEList object

# Filtering -----------------------------------------------------------------------------------
index.keep.expr <- filterByExpr(dgelist, group = sample.meta.data$Tissue) # FIXED: Correct grouping

# Validate filtering to ensure appropriate number of genes retained
print(table(index.keep.expr)) # FIXED: Added for validation

# Filter lowly expressed genes ----------------------------------------------------------------
dgelist.filtered <- dgelist[index.keep.expr, , keep.lib.sizes = FALSE]

# Normalization and filtering on filtered list -----------------------------------------------
dgelist.filtered.norm <- calcNormFactors(dgelist.filtered, method = "TMM")

# Calculate LCPM of new filtered list --------------------------------------------------------
lcpm.filtered <- cpm(dgelist.filtered, log = TRUE)

# Create a long (tidy) format data frame for LogCPM values ------------------------------------
long.lcpm.filtered <- as.data.frame(lcpm.filtered) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogPCM", names_to = "SRR_ID")

# Create a density plot of the sample distribution -------------------------------------------
ggplot(long.lcpm.filtered) + geom_density(aes(LogPCM, color = SRR_ID)) +
  labs(title = "Log CPM Distribution for Filtered Genes",
       x = "Log CPM",
       y = "Density")

# Combining CPM and LogCPM into a single dataframe ===========================================
# Calculate CPM and Log2 CPM values from the dgelist.filtered.norm object --------------------
cpm.dgelist.filtered.norm <- cpm(dgelist.filtered.norm)
lcpm.dgelist.filtered.norm <- cpm(dgelist.filtered.norm, log = TRUE)

# Boxplots to show before and after normalization ============================================
# Calculate LogCPM for raw counts (before normalization, use filtered data) ------------------
raw.lcpm <- cpm(dgelist.filtered, log = TRUE) # FIXED: Use filtered data for consistency

# Create a long format dataframe for raw LogCPM values (before normalization)
long.raw.lcpm <- as.data.frame(raw.lcpm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID") %>%
  mutate(Status = "Before Normalization")

# Create a long format dataframe for normalized LogCPM values (after normalization)
long.norm.lcpm <- as.data.frame(lcpm.dgelist.filtered.norm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID") %>%
  mutate(Status = "After Normalization")

# Combine both data frames for plotting -----------------------------------------------------
combined.lcpm <- bind_rows(long.raw.lcpm, long.norm.lcpm)

# Plot the boxplots --------------------------------------------------------------------------
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
