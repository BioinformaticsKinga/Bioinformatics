import pandas as pd
from matplotlib_venn import venn3_unweighted
import matplotlib.pyplot as plt

# Step 1: Load gene expression results from three different methods
res_deseq2 = pd.read_csv('deseq2_results.txt', sep='\t')  # Read DESeq2 results
res_edger = pd.read_csv('edger_results.txt', sep='\t')    # Read edgeR results
res_limma = pd.read_csv('limma_results.txt', sep='\t')    # Read limma-voom results

# Step 2: Extract unique gene sets from each method
genes_deseq2 = set(res_deseq2['GeneID'])  # Convert DESeq2 gene list to a set
genes_edger = set(res_edger['GeneID'])    # Convert edgeR gene list to a set
genes_limma = set(res_limma['GeneID'])    # Convert limma-voom gene list to a set

# Step 3: Compute the intersection of all three sets (common genes)
common_genes = genes_deseq2 & genes_edger & genes_limma  # Find common genes
num_common_genes = len(common_genes)  # Count the number of common genes
print(f"Number of common genes: {num_common_genes}")  # Print the result

# Step 4: Create a Venn diagram (equal-sized circles)
plt.figure(figsize=(8, 8))  # Set figure size
venn_diagram = venn3_unweighted(
    [genes_deseq2, genes_edger, genes_limma], 
    ('DESeq2', 'edgeR', 'limma_voom')  # Labels for the three methods
)

# Step 5: Set colors for different sections of the Venn diagram
venn_diagram.get_patch_by_id('100').set_color('tomato')       # Red for DESeq2 only
venn_diagram.get_patch_by_id('010').set_color('royalblue')    # Blue for edgeR only
venn_diagram.get_patch_by_id('001').set_color('aqua')         # Aqua for limma-voom only

# Set colors for overlapping regions
venn_diagram.get_patch_by_id('110').set_color('lightcoral')   # Light red for DESeq2 + edgeR
venn_diagram.get_patch_by_id('101').set_color('lightcoral')   # Light red for DESeq2 + limma
venn_diagram.get_patch_by_id('011').set_color('skyblue')      # Light blue for edgeR + limma
venn_diagram.get_patch_by_id('111').set_color('lightgray')    # Gray for all three methods

# Step 6: Add black borders around the circles
for patch in venn_diagram.patches:
    patch.set_edgecolor('black')  # Set border color to black
    patch.set_linewidth(5)        # Increase border thickness

# Step 7: Adjust font size and boldness of numbers inside the circles
for text in venn_diagram.subset_labels:  # Labels inside the circles (numbers)
    if text:  # Ensure text is not None
        text.set_fontweight('bold')  # Make numbers bold
        text.set_fontsize(16)        # Increase font size

# Step 8: Adjust font size and boldness of category labels (DESeq2, edgeR, limma)
for text in venn_diagram.set_labels:  # Labels for the three sets
    if text:  # Ensure text is not None
        text.set_fontweight('bold')  # Make labels bold
        text.set_fontsize(18)        # Increase font size

# Step 9: Save and display the final Venn diagram
plt.savefig("venn_diagram_epip24.png")  # Save as PNG file
plt.show()  # Display the diagram
