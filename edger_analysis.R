# Load necessary libraries
library(edgeR)

# Load count data and experimental conditions
counts <- read.table("counts.txt", header=TRUE, row.names=1)
coldata <- data.frame(condition=c("control", "treatment"))  # Replace with actual conditions

# Prepare the data for edgeR
y <- DGEList(counts=counts, group=coldata$condition)
y <- calcNormFactors(y)

# Fit the model and perform the statistical test
design <- model.matrix(~ coldata$condition)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Perform the test
qlf <- glmQLFTest(fit)

# Obtain results
res <- topTags(qlf, n=Inf)

# Write results to CSV
write.csv(as.data.frame(res), file="zygmukin_edger_results.csv")
