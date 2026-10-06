# ============================================================

# RNA-seq differential expression analysis

# Nitisinone treatment in Aedes aegypti

# ============================================================

#

# Purpose:

# Differential expression analysis of RNA-seq data generated

# from Aedes aegypti following nitisinone administration.

#

# Experimental design:

# A 2 x 2 factorial design was used with:

#

# condition: Control vs Nitisinone

# time:      24 h vs 48 h

#

# Nitisinone was administered together with a blood meal.

# Control insects received a blood meal containing PBS.

#

# The experimental groups consisted of independent insects.

#

# Main analyses:

#

# 1. Filtering of low-expression/high-variability genes

# using HTSFilter.

#

# 2. rlog transformation for exploratory analyses.

#

# 3. Differential expression analysis using DESeq2 with a

# full factorial model:

#

# ~ condition + time + condition:time

#

# 4. Likelihood Ratio Test (LRT) to test the

# condition-by-time interaction.

#

# 5. Treatment-specific comparisons at 24 h and 48 h.

#

# 6. Log2 fold-change shrinkage using apeglm.

#

# ============================================================

# ============================================================

# 1. Software requirements

# ============================================================

# Required R packages:

#

# DESeq2

# HTSFilter

# apeglm

#

# Install packages if necessary:

#

# if (!requireNamespace("BiocManager", quietly = TRUE))

# install.packages("BiocManager")

#

# BiocManager::install(c(

# "DESeq2",

# "HTSFilter",

# "apeglm"

# ))

# ============================================================

# 2. Load packages

# ============================================================

library(DESeq2)
library(HTSFilter)
library(apeglm)

# ============================================================

# 3. Define input and output directories

# ============================================================

# The repository should contain separate directories for input

# data and analysis results.

#

# Update these paths according to the location of the repository

# on your system.

input_dir <- "input"
output_dir <- "results"

# Create the output directory if it does not already exist.

if (!dir.exists(output_dir)) {
dir.create(output_dir, recursive = TRUE)
}

# ============================================================

# 4. Load raw count matrix

# ============================================================

# The count matrix must contain:

#

# - one column with gene identifiers

# - one column per RNA-seq sample

#

# The sample names in the count matrix must match the sample

# identifiers in the experimental metadata file.

counts_file <- file.path(
input_dir,
"Counts_all.txt"
)

Counts_all <- read.delim(
counts_file,
sep = "\t",
header = TRUE,
stringsAsFactors = FALSE,
check.names = FALSE
)

# Extract gene identifiers.

gene_ids <- Counts_all[, 1]

# Convert the count table to a numeric matrix.

Counts_all_matrix <- data.matrix(
Counts_all[, -1]
)

rownames(Counts_all_matrix) <- gene_ids

# ============================================================

# 5. Load experimental metadata

# ============================================================

# The metadata file must contain at least:

#

# sample

# condition

# time

#

# condition:

# Control   = PBS-containing blood meal

# Treated   = nitisinone-containing blood meal

#

# time:

# 24        = 24 h after treatment

# 48        = 48 h after treatment

metadata_file <- file.path(
input_dir,
"Experiment_description_all.txt"
)

experiment_all <- read.delim(
metadata_file,
sep = "\t",
header = TRUE,
stringsAsFactors = FALSE,
check.names = FALSE
)

# ============================================================

# 6. Validate count matrix and metadata

# ============================================================

# Check that all samples in the metadata are present in the

# count matrix.

if (!all(experiment_all$sample %in% colnames(Counts_all_matrix))) {
missing_samples <- setdiff(
experiment_all$sample,
colnames(Counts_all_matrix)
)

```
stop(
    "The following samples are present in the metadata ",
    "but missing from the count matrix: ",
    paste(missing_samples, collapse = ", ")
)
```

}

# Check that the number of samples is identical.

if (ncol(Counts_all_matrix) != nrow(experiment_all)) {
stop(
"The number of samples in the count matrix does not ",
"match the number of rows in the metadata."
)
}

# Reorder the count matrix to match the metadata.

Counts_all_matrix <- Counts_all_matrix[
,
experiment_all$sample,
drop = FALSE
]

# Confirm that sample names are now identical and in the same order.

stopifnot(
identical(
colnames(Counts_all_matrix),
experiment_all$sample
)
)

# ============================================================

# 7. Define experimental factors

# ============================================================

# Set factor levels explicitly to define the reference groups.

#

# Control is the reference level for condition.

# 24 h is the reference level for time.

experiment_all$condition <- factor(
experiment_all$condition,
levels = c("Control", "Treated")
)

experiment_all$time <- factor(
experiment_all$time,
levels = c("24", "48")
)

# Check the experimental design.

print(
table(
experiment_all$condition,
experiment_all$time
)
)

# ============================================================

# 8. Construct the DESeq2 dataset

# ============================================================

# Full factorial model:

#

# expression ~ condition + time + condition:time

#

# This model estimates:

#

# - the main effect of treatment

# - the main effect of time

# - the treatment-by-time interaction

dds_all <- DESeqDataSetFromMatrix(
countData = Counts_all_matrix,
colData = experiment_all,
design = ~ condition + time + condition:time
)

# ============================================================

# 9. Filter low-expression / highly variable genes

# ============================================================

# HTSFilter is used to remove genes with low expression and/or

# excessive variability across samples.

#

# s.len = 100 corresponds to the number of candidate models

# evaluated by HTSFilter.

dds_filtered <- HTSFilter(
dds_all,
s.len = 100,
plot = TRUE
)$filteredData

# Save the number of genes retained after filtering.

cat(
"Genes before filtering:",
nrow(dds_all),
"\n"
)

cat(
"Genes after HTSFilter:",
nrow(dds_filtered),
"\n"
)

# ============================================================

# 10. rlog transformation

# ============================================================

# rlog-transformed expression values are generated for

# exploratory analyses and visualization.

#

# The transformed data are NOT used for DESeq2 differential

# expression testing.

rlog_all <- rlog(
dds_filtered,
blind = TRUE,
fitType = "parametric"
)

# Extract transformed expression values.

rlog_matrix <- assay(
rlog_all
)

# Save rlog-transformed expression matrix.

write.table(
rlog_matrix,
file = file.path(
output_dir,
"RLOG_transformed_values_all.txt"
),
sep = "\t",
dec = ".",
row.names = TRUE,
col.names = NA,
quote = FALSE
)

# ============================================================

# 11. Fit the full DESeq2 model

# ============================================================

dds_complete <- DESeq(
dds_filtered
)

# Display the coefficients estimated by DESeq2.

coef_names <- resultsNames(
dds_complete
)

print(coef_names)

# ============================================================

# 12. Likelihood Ratio Test (LRT)

# ============================================================

# The LRT tests the condition-by-time interaction.

#

# Full model:

#

# ~ condition + time + condition:time

#

# Reduced model:

#

# ~ condition + time

#

# Therefore, the LRT tests whether the treatment effect differs

# between 24 h and 48 h.

dds_LRT <- DESeq(
dds_complete,
test = "LRT",
reduced = ~ time + condition
)

# Extract LRT results.

results_LRT <- results(
dds_LRT,
independentFiltering = FALSE
)

# Save the complete LRT results.

write.table(
as.data.frame(results_LRT),
file = file.path(
output_dir,
"DESeq2_LRT_condition_by_time_interaction.txt"
),
sep = "\t",
dec = ".",
row.names = TRUE,
col.names = NA,
quote = FALSE
)

# ============================================================

# 13. Treatment effect at 24 h

# ============================================================

# Subset the DESeq2 dataset to samples collected 24 h after

# treatment.

dds_24h <- dds_complete[
,
colData(dds_complete)$time == "24"
]

# For this subset, the analysis evaluates only the treatment

# effect.

design(dds_24h) <- ~ condition

# Refit the model using only the 24 h samples.

dds_24h <- DESeq(
dds_24h
)

# Define the specific comparison:

#

# Treated (nitisinone) vs Control (PBS)

contrast_24h <- c(
"condition",
"Treated",
"Control"
)

# Obtain differential expression results.

res_24h <- results(
dds_24h,
contrast = contrast_24h,
independentFiltering = FALSE
)

# Shrink log2 fold changes using apeglm.

#

# lfcThreshold = 1 specifies a log2 fold-change threshold of 1.

resLFC_24h <- lfcShrink(
dds_24h,
coef = "condition_Treated_vs_Control",
res = res_24h,
type = "apeglm",
lfcThreshold = 1
)

# Save 24 h results.

write.table(
as.data.frame(resLFC_24h),
file = file.path(
output_dir,
"DESeq2_statistics_24h_Treated_vs_Control.txt"
),
sep = "\t",
dec = ".",
row.names = TRUE,
col.names = NA,
quote = FALSE
)

# ============================================================

# 14. Treatment effect at 48 h

# ============================================================

# Subset the DESeq2 dataset to samples collected 48 h after

# treatment.

dds_48h <- dds_complete[
,
colData(dds_complete)$time == "48"
]

# For this subset, the analysis evaluates only the treatment

# effect.

design(dds_48h) <- ~ condition

# Refit the model using only the 48 h samples.

dds_48h <- DESeq(
dds_48h
)

# Define the specific comparison:

#

# Treated (nitisinone) vs Control (PBS)

contrast_48h <- c(
"condition",
"Treated",
"Control"
)

# Obtain differential expression results.

res_48h <- results(
dds_48h,
contrast = contrast_48h,
independentFiltering = FALSE
)

# Shrink log2 fold changes using apeglm.

resLFC_48h <- lfcShrink(
dds_48h,
coef = "condition_Treated_vs_Control",
res = res_48h,
type = "apeglm",
lfcThreshold = 1
)

# Save 48 h results.

write.table(
as.data.frame(resLFC_48h),
file = file.path(
output_dir,
"DESeq2_statistics_48h_Treated_vs_Control.txt"
),
sep = "\t",
dec = ".",
row.names = TRUE,
col.names = NA,
quote = FALSE
)

# ============================================================

# 15. End of analysis

# ============================================================

cat("\n")
cat("RNA-seq differential expression analysis completed.\n")
cat(
"Results were written to:",
normalizePath(output_dir),
"\n"
)

