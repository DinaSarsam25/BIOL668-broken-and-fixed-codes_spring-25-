

```{r}
library( "DESeq2" )
library(tidyverse)
library(ggplot2)
library(dplyr)


```

```{r}
count_data_file <- 'data/rna_counts_data.csv'
map_file <- 'data/rna_map_updatecopy.csv'
# the first two codes are showing the files paths.
```
```{r}
counts_path <- "C:/Users/MyPC/Downloads/rna_counts_data.csv"
map_path <- "C:/Users/MyPC/Downloads/rna_map_update copy.csv"

#the two above codes are setting a new paths for the computer.


rna_counts <- read.csv(counts_path, header = TRUE, row.names = 1, fileEncoding = "UTF-8-BOM")
rna_map <- read.csv(map_path, header = TRUE, fileEncoding = "UTF-8-BOM")
# read.csv code used to download the file
#header = true to contains column names
#row.names = 1 --> means that the first column in the CSV files is the first row name. Which is the gene name.
#fileEncoding = "UTF-8-BOM" --> read the order correctly. 

head(rna_counts)
head(rna_map)
```
#the rna_counts data is showing: first column is showing gene names while the second is showing the length of gene in base pairs. 
#from the third column to 23 are represent the name of a sample and the number under this column is the raw read counts for each gene in the sample.

#rna map is showing the first column that mean sample.number which is ID used. Cloumn two showing sample.group 

#tumor.type and mutation group.

```{r}
setdiff(rna_map$Sample, colnames(rna_counts)) #data processing step.
rna_map <- rna_map[match(colnames(rna_counts), rna_map$Sample), ]
all(rna_map$Sample == colnames(rna_counts))

# Check the first few row names of rna_counts to make sure they represent genes
head(rownames(rna_counts))

# Check if any row names are NA or incorrect
any(is.na(rownames(rna_counts)))
rna_counts <- read.csv(counts_path, header = TRUE, row.names = 1, fileEncoding = "UTF-8-BOM")
```
```{r}
# Load necessary libraries
library(DESeq2)

# Set file paths for count data and metadata and code for reads a CSV file into R.
counts_path <- "C:/Users/MyPC/Downloads/rna_counts_data.csv"
map_path <- "C:/Users/MyPC/Downloads/rna_map_update copy.csv"

# To Read in the count data and metadata
rna_counts <- read.csv(counts_path, header = TRUE, row.names = 1, fileEncoding = "UTF-8-BOM") #This tells R that the first column contains row names.
rna_map <- read.csv(map_path, header = TRUE, fileEncoding = "UTF-8-BOM")  #UTF-8-BOM": Specifies the file encoding

# To Check the first few rows of rna_counts and rna_map to inspect the data
head(rna_counts)
head(rna_map)

# To Ensure column names in rna_counts match sample names in rna_map
# Check column names in rna_counts and the sample identifiers in rna_map
colnames(rna_counts)
head(rna_map$Sample.Number)

# To reorganize rna_map to match the column names of rna_counts (sample identifiers)
rna_map <- rna_map[match(colnames(rna_counts), rna_map$Sample.Number), ]

# Check that the reorder was successful
all(rna_map$Sample.Number == colnames(rna_counts))  # 
```
```{r}
# To Ensure the Tumor.Type column is a factor and inspect the levels
rna_map$Tumor.Type <- factor(rna_map$Tumor.Type)
table(rna_map$Tumor.Type)  # Check levels like R132H, R132Q, WT
rna_map <- as.data.frame(rna_map)  # To Convert to data frame 
dim(rna_counts)
dim(rna_map)
rownames(rna_map) <- colnames(rna_counts)
rownames(rna_map) <- colnames(rna_counts)
colnames(rna_map)
if (!all(rownames(rna_map) %in% colnames(rna_counts))) {
  stop("Row names of rna_map must match column names of rna_counts.")   #(I added this code because i do not know why sudnly my code was not working, but i couldnot figure it out)

}
dim(rna_counts)
dim(rna_map)
rownames(rna_map) <- colnames(rna_counts)
rownames(rna_map) <- colnames(rna_counts)
colnames(rna_map)
# To Create DESeqDataSet from the count matrix and metadata

dds <- DESeqDataSetFromMatrix(countData = rna_counts, 
                              rna_counts_1 = sample.number
                              rna_map = sample_data



# results for the comparison (e.g., R132H vs WT)
res <- results(deSeqDataSet, contrast = c("Tumor.Type", "R132H", "WT"))

# View the results
head(res)

# Save the results to a CSV file
write.csv(as.data.frame(res), "differential_expression_results.csv")

# Plot MA plot to visualize differential expression
plotMA(res, ylim = c(-5, 5))

# Optional: Volcano plot for visualizing differential expression

volcano_data <- as.data.frame(res)
volcano_data$significant <- ifelse(volcano_data$padj < 0.05, "Yes", "No")
plot(volcano_data$log2FoldChange, -log10(volcano_data$padj), 
     pch = 20, col = ifelse(volcano_data$significant == "Yes", "red", "black"), 
     xlab = "Log2 Fold Change", ylab = "-Log10 p-value")
abline(h = -log10(0.05), col = "blue", lty = 2)  # Add significance threshold line

```
