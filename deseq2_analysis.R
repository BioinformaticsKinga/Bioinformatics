# Load necessary libraries
library(DESeq2)

# Load count data and sample information
counts <- read.table("counts.txt", header=TRUE, row.names=1)
coldata <- data.frame(condition=c("control", "treatment"))  # Replace with actual sample conditions

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# Perform DESeq2 analysis
dds <- DESeq(dds)

# Obtain results
res <- results(dds)

# Write results to a CSV file
write.csv(as.data.frame(res), file="zygmukin_deseq2_results.csv")
