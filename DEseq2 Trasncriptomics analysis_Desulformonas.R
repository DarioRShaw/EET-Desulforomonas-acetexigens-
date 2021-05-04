
#Transcriptomics analysis Deseq2 of Desulforomonas project
#After have the count table (Genes vs raw counts)

#Install DEseq2
source("https://bioconductor.org/biocLite.R")
biocLite()  
biocLite("DESeq2")

#Load the library
library('DESeq2') 

#Other libraries
library('RColorBrewer')
library (ggplot2) 
library(dplyr)
library("gplots")



#let's see a general overview of the data 
#General information of the data
#inport data 
countData = read.table(file = "~/Desktop/Dario Rangel Shaw/ANAMMOX research/Desulforomonas project/Bioinformatic analysis/Bioinformatic FINAL analysis for manuscript /DEG analysis/Input files for DESeq2/Count_table_fumarate_vs_set_potential_vs_OCV.txt", header = TRUE, sep = '\t', row.names = 1) 
dim(countData)
head(countData)

#Create experiments conditions labels
colData = read.table(file = "~/Desktop/Dario Rangel Shaw/ANAMMOX research/Desulforomonas project/Bioinformatic analysis/Bioinformatic FINAL analysis for manuscript /DEG analysis/Input files for DESeq2/colData_fumarate_vs_set_potential_vs_OCV.txt", header = TRUE, sep = '\t', row.names = 1) 

#Run DEseq algorith
dds <- DESeqDataSetFromMatrix(countData, colData = colData, design = ~ condition)
dds <- DESeq(dds)
dds <- DESeq(dds, minReplicatesForReplace=Inf)
res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)


#First let's visualize with a PCA the different evaluated conditions and replicates
#PCA analysis of the samples (group according expression profile between samples)
#The differential expression analysis above operates on the raw (normalized) count data. 
#But for visualizing or clustering data  you ned to work with transformed versions of the data. 
#First, use a regularlized log trasnformation while re-estimating the dispersion ignoring any information
#you have about the samples (blind=TRUE). Perform a principal components analysis and hierarchical clustering.
#rld <- rlogTransformation(dds, blind=TRUE)
#plotPCA(rld, intgroup= 'condition' ) 

#PCA with ggplot2
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1,PC2, color=condition,shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


#Now let's see a Hierarchical clustering analysis 
#let's get the actual values for the first few genes 
head(assay(rld))
## now transpose those
t(head(assay(rld)))
## now get the sample distances from the transpose of the whole thing
dist(t(assay(rld)))
sampledist <- dist(t(assay(rld)))
plot(hclust(sampledist))


#We can also visualize a Heatmap with the hierachical clustering
sampledist
as.matrix(sampledist)
sampledistmat <- as.matrix(sampledist)
heatmap(sampledistmat)

#Let's make it more fancy
Colors=colorRampPalette(brewer.pal(10, "RdBu"))(300)
heatmap.2(sampledistmat, density.info = "none",trace = "none",col=Colors, sepwidth=c(0.0005, 0.0005), sepcolor ="black",colsep=1:nrow(sampledistmat),rowsep=1:nrow(sampledistmat)) 

#with viridis palette
heatmap.2(sampledistmat, density.info = "none",trace = "none",col=viridis::viridis_pal(), sepwidth=c(0.0005, 0.0005), sepcolor ="black", colsep=1:nrow(sampledistmat),rowsep=1:nrow(sampledistmat)) 

#What about Spectral palette
#pdf(file = "Heat map with hierarchical clustering.pdf")
Colors=colorRampPalette(brewer.pal(10, "Spectral"))(300)
heatmap.2(sampledistmat, density.info = "none",trace = "none",col=Colors, sepwidth=c(0.00005, 0.00005), sepcolor ="black", colsep=1:nrow(sampledistmat),rowsep=1:nrow(sampledistmat)) 
#dev.off()
#I liked more this one, so let's use this palette to save it



#Now let's perform the DESeq2 analysis to evaluate the difference between Set potential vs Fumarate vs OCV 
#First we will evaluate between the DEG from Fumarate to Set potential
#For DESeq2 analysis first we import the data
#Import the count table
countData = read.table(file = "~/Desktop/Dario Rangel Shaw/ANAMMOX research/Desulforomonas project/Bioinformatic analysis/Bioinformatic FINAL analysis for manuscript /DEG analysis/Input files for DESeq2/Count_table_fumarate_vs_set_potential.txt", header = TRUE, sep = '\t', row.names = 1) 

#Check that Cout Table was imported properly
dim(countData)
head(countData)

#Create experiments conditions labels
colData = read.table(file = "~/Desktop/Dario Rangel Shaw/ANAMMOX research/Desulforomonas project/Bioinformatic analysis/Bioinformatic FINAL analysis for manuscript /DEG analysis/Input files for DESeq2/colData_fumarate_vs_set_potential.txt", header = TRUE, sep = '\t', row.names = 1) 

#Create DEseq input matrix
dds <- DESeqDataSetFromMatrix(countData, colData = colData, design = ~ condition)

#Run DEseq algorith
dds <- DESeq(dds)

#Sometimes you can get Padj values of 0, to know more about the reasons, you can read 
#https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
#To have the only p values set to NA with those from genes with all counts equal to zero we execute
dds <- DESeq(dds, minReplicatesForReplace=Inf)

# Get differential expression results
#res <- results(dds)
res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
head(res)

#let's see and verify if our P-values has a proper distribution
#Histogram p-values
hist(res$pvalue, breaks=50, col="grey")
#Everything is good

#Now let's check if our genes fits our model
#Dispersion plot 
plotDispEsts(dds)
#Everything is good

#visualize differentially expressed genes in MA plots
#svg(file = "MA plot Fumarate vs Set potential.svg")
plotMA(dds)
#dev.off()

#Optional
#order by adjusted p-value
#res <- res[order(res$padj), ]
#head(res)

#Export results of differential analysis expression in a table 
write.csv(res, file = "Desulforomonas DESEq2 analysis Fumarate vs Set potential.csv") 

#Optional
#Export just genes with significant change in expression Padj <0.05
#sig <- subset(res, padj < 0.05)
#write.csv(sig, file = "Desulforomonas DESEq2 analysis Fumarate vs Set potential.csv")  



#Now we will evaluate between the DEG from OCV to Set potential
#For DESeq2 analysis first we import the data
#Import the count table
countData = read.table(file = "~/Desktop/Dario Rangel Shaw/ANAMMOX research/Desulforomonas project/Bioinformatic analysis/Bioinformatic FINAL analysis for manuscript /DEG analysis/Input files for DESeq2/Count_table_OCV_vs_set_potential.txt", header = TRUE, sep = '\t', row.names = 1) 

#Check that Cout Table was imported properly
dim(countData)
head(countData)

#Create experiments conditions labels
colData = read.table(file = "~/Desktop/Dario Rangel Shaw/ANAMMOX research/Desulforomonas project/Bioinformatic analysis/Bioinformatic FINAL analysis for manuscript /DEG analysis/Input files for DESeq2/colData_OCV_vs_set_potential.txt", header = TRUE, sep = '\t', row.names = 1) 

#Create DEseq input matrix
dds <- DESeqDataSetFromMatrix(countData, colData = colData, design = ~ condition)

#Run DEseq algorith
dds <- DESeq(dds)

#Sometimes you can get Padj values of 0, to know more about the reasons, you can read 
#https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
#To have the only p values set to NA with those from genes with all counts equal to zero we execute
dds <- DESeq(dds, minReplicatesForReplace=Inf)

# Get differential expression results
#res <- results(dds)
res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
head(res)

#let's see and verify if our P-values has a proper distribution
#Histogram p-values
hist(res$pvalue, breaks=50, col="grey")
#Everything is good

#Now let's check if our genes fits our model
#Dispersion plot 
plotDispEsts(dds)
#Everything is good

#visualize differentially expressed genes in MA plots
#svg(file = "MA plot OCV vs Set potential.svg")
plotMA(dds)
#dev.off()

#Optional
#order by adjusted p-value
#res <- res[order(res$padj), ]
#head(res)

#Export results of differential analysis expression in a table 
write.csv(res, file = "Desulforomonas DESEq2 analysis OCV vs Set potential.csv") 

#Optional
#Export just genes with significant change in expression Padj <0.05
#sig <- subset(res, padj < 0.05)
#write.csv(sig, file = "Desulforomonas DESEq2 analysis OCV vs Set potential.csv")  







  
