Part 1: Tutorial using limma, Glimma and Edge R using the given datasets.
Project goal is to analyse RNA-sequencing data from the mouse mammary gland using three packages:
1-edgeR--> import, organise, filter and normalise the dat 
2-Limma --> To assess differential expression and perform gene set testing. It consist of linear modelling and empirical Bayes moderation. 
3-Glimma --> To enables interactive exploration of the results so that individual samples and genes can be examined by the user.

 
```{r}
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(AnnotationDbi)

# Install the 'airway' package from Bioconductor ==> BiocManager is installed to manage Bioconductor packages.The airway package is installed from Bioconductor, which is necessary for performing RNA-seq analysis and is commonly used in bioinformatics tutorials.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("airway")
```
```{r}

url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
   "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
   "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
   "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
   "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)

# The code above used to Downloads the GSE63310 dataset from NCBI GEO as a tarball file. Extracts the contents of the tarball into the current working directory. Unzips each of the compressed .gz files that are part of the dataset. Then, Reads and inspects the first few rows of the first data file to verify the format and contents.
```

```{r}
x <- readDGE(files, columns=c(1,3)) # to read the data
class(x)  # => a specific class used in edgeR to store RNA-seq count data and associated information.
dim(x) # this line to returns the dimensions of the object 
```
Organising sample information
```{r}
samplenames <- substring(colnames(x), 12, nchar(colnames(x))) # Cleaning up sample names
samplenames
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP")) #this code to access the experimental group.
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2))) #This creates a factor called lane to represent which sequencing lane each sample came from.
x$samples$lane <- lane 
x$samples
```
Organising gene annotations
```{r}
geneid <- rownames(x) # grabs the gene IDs from the DGEList object x
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)
genes <- genes[!duplicated(genes$ENTREZID),] # to remove Remove Duplicate Gene Annotations
x$genes <- genes
x
```
Data pre-processing
5.1Transformations from the raw-scale

```{r}
cpm <- cpm(x) #cpm() function from the edgeR package to calculate normalized expression values for each gene in each sample. 
lcpm <- cpm(x, log=TRUE) #Visualizations (like heatmaps, PCA)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
```
5.2Removing genes that are lowly expressed
```{r}
table(rowSums(x$counts==0)==9)
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
```
```{r}
lcpm.cutoff <- log2(10/M + 2/L) #calculate a Suggested log-CPM Cutoff
library(RColorBrewer) #Set Colors for Each Sample
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2)) #Set up Side-by-Side Plots
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="") #Plot Density of log-CPM for All Samples
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}             #
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")  #re-calculate this after filtering lowly expressed genes
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
```
5.3Normalising gene expression distributions
```{r}
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors   # to ensure fair comparisons between samples by accounting for technical biases.
```
```{r}
x2 <- x
x2$samples$norm.factors <- 1  # to duplicate the Original Dataset
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5  #Simulates severe differences in library composition
par(mfrow=c(1,2))          #Calculates log-CPM without normalization and boxplot for each sample. 
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)      #og-CPM values reflect normalized expression
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
```
5.4Unsupervised clustering of samples
```{r}
lcpm <- cpm(x, log=TRUE)   #Compute log-CPM
par(mfrow=c(1,2))
col.group <- group  #Color Samples by Group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane        #Color Samples by Sequencing Lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)  #Create MDS Plot for Sample Groups
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)  #3D-like Interactive MDS Plot 
```
6Differential expression analysis
6.1Creating a design matrix and contrasts
```{r}
design <- model.matrix(~0+group+lane) #Builds a design matrix adjusting for biological group and sequencing batch.
colnames(design) <- gsub("group", "", colnames(design))
design
contr.matrix <- makeContrasts(       #Defines pairwise comparisons of interest between biological conditions
   BasalvsLP = Basal-LP, 
   BasalvsML = Basal - ML, 
   LPvsML = LP - ML, 
   levels = colnames(design))
contr.matrix
```
6.2Removing heteroscedascity from count data
```{r}
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE) #Normalize and estimate variance per gene
v
vfit <- lmFit(v, design) #Fit a linear model per gene
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) #Apply comparisons (e.g. Basal vs LP)
efit <- eBayes(vfit) #Statistical testing with empirical Bayes
plotSA(efit, main="Final model: Mean-variance trend") #Visualize mean-variance trend of final model.
```
6.4Examining the number of DE genes
```{r}
summary(decideTests(efit))
tfit <- treat(vfit, lfc=1) #filters DE genes based on meaningful fold-change
dt <- decideTests(tfit) #Gets significant DE results per contrast
summary(dt)
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20) #Summarizes DE results per group
write.fit(tfit, dt, file="results.txt") #Exports DE results to file
```
6.5Examining individual DE genes from top to bottom

```{r}
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf) #op differentially expressed genes based on the treat model.
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf) #extracting results from the first contrast
head(basal.vs.lp)
head(basal.vs.ml)
```
6.6Useful graphical representations of differential expression results
```{r}
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=FALSE)  #To create a Mean-Difference plot to visualizing the relationship between the log fold change (x-axis) and mean expression (y-axis) of genes in the contrast.

library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")


heatmap.2(lcpm[i,], scale="row",   #from the gplots package is used to create a heatmap of the top 100 differentially expressed genes.
   labRow=v$genes$SYMBOL[i], labCol=group, 
   col=mycol, trace="none", density.info="none", 
   margins =c(5,5),cexRow = 0.7,           
   cexCol = 0.8,lhei=c(1,4), dendrogram="column")
```
7Gene set testing with camera

```{r}


load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123")) #pre-defined mouse gene set


idx <- ids2indices(Mm.c2,id=rownames(v))
cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.BasalvsLP,5) #pathways or gene sets used for pathway enrichment analysis.
cam.BasalvsML <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.BasalvsML,5) #Camera for Basal vs LP
cam.LPvsML <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.LPvsML,5)
barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, 
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")
```
8Software and code used
```{r}
sessionInfo() #to get a summary about the codes in this notebook.
```
