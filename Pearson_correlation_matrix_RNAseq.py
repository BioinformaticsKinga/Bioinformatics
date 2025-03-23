import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load RNA-Seq data
# Assuming that the data has genes in rows and samples in columns
data = pd.read_csv("rna_seq_data_lupinus.csv", index_col=0)

# Compute Pearson correlation between the samples (columns)
correlation_matrix = data.corr(method="pearson")

# Round the correlation values to 3 decimal places
correlation_matrix = correlation_matrix.round(3)

# Print the correlation matrix
print(correlation_matrix)

# Create a heatmap to visualize the Pearson correlation between the samples
plt.figure(figsize=(10, 8))
sns.set(font_scale=0.8)  # Adjusting font size
sns.heatmap(correlation_matrix, annot=True, fmt=".3f", cmap="viridis", 
            linewidths=0.5, cbar=True, center=0, square=True)
plt.title("Pearson Correlation Heatmap Between Samples")
plt.show()
