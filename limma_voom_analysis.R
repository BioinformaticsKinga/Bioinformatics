# Load necessary libraries
library(limma)
library(edgeR)

# Load count data
counts <- read.table("counts.txt", header=TRUE, row.names=1)
coldata <- data.frame(condition=c("control", "treatment"))  

# Prepare data for voom transformation
d <- DGEList(counts=counts, group=coldata$condition)
d <- calcNormFactors(d)
v <- voom(d, design=model.matrix(~ coldata$condition), plot=TRUE)

# Fit the model and perform eBayes
fit <- lmFit(v, design=model.matrix(~ coldata$condition))
fit <- eBayes(fit)

# Obtain results
res <- topTable(fit, number=Inf)

# Write results to CSV file
write.csv(res, file="zygmukin_limma_voom_results.csv")
