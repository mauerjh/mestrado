---
title: "Quality Control of miRNA sequencing (Code/Complete Version)"
author: "Jessica"
date: "10/09/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

```

We will need two files: 1) sample information table (`coldata`), and 2) count matrix data (`cts`) And format them so data is in the same order.

```{r message=FALSE, warning=FALSE}
library("DESeq2")
library ("ggplot2")
library("ggfortify")
library('EnhancedVolcano')
library("ggrepel")
library("dplyr")

cts = read.table(file = "Samples-maturemiRNAcounts.tsv", header = T)
dim(cts)
coldata = read.table(file = "info_sex_age.txt", header = T)
coldata$condition <- rep(c("control", "case"), each=12)
dim(coldata)

# Retirada da amostra que não sequenciou e das amostras teste
unicount <- read.csv("4b_ReadCountUniquely.txt", row.names=1, sep="")
uni_inpd <- select (unicount, -8, -25, -26, -27, -28)

#Preparacao da tabela de readcounts (countData)
mRC_uni <- as.matrix(uni_inpd)
mRC_uni <- as.data.frame (mRC_uni)
colnames(mRC_uni) <- gsub("uni", "", colnames(mRC_uni))
colnames(mRC_uni) <- gsub("X", "", colnames(mRC_uni))

#reading sample information file 
sampleinfo <- read.csv2("sample_information_140521.csv", h=T)
#checar e arrumar o que e factor, numeric etc
#deixar o id igual ao da count data
sampleinfo$ID <- gsub("_", ".", sampleinfo$ID)

#deixar a sampleinfo na mesma ordem que o count data
as.data.frame(colnames(mRC_uni)) -> ids
colnames(ids) <- "ID"
join(ids, sampleinfo, by="ID", type="inner") -> sampleinfo_order
rownames(sampleinfo_order) <- sampleinfo_order$ID

#checar se a ordem esta igual
all(rownames(sampleinfo_order) == colnames(mRC_uni))

#tirar a primeira coluna de id
sampleinfo_todosgroup <- sampleinfo_order[, c(-2,-3)]

```


First, we can look at total read counts per sample
```{r include = TRUE, message=FALSE, warning=FALSE}
barplot(colSums(mRC_uni), las=3,main="Total uniquely mapped miRNA read counts per sample", ylab="Total miRNA read counts", border = NA)

```

Then, we can look to mean read counts per gene. First we can filter only genes with higher expression (the first 30)
```{r include = TRUE, message=FALSE, warning=FALSE}
cts_order <- mRC_uni[order(rowMeans(mRC_uni), decreasing = T),]
cts_order_f <- cts_order[1:30,]
barplot(rowMeans(cts_order_f[,]), las=2,main="Mean read counts per gene", ylab="Mean read counts",
         cex.axis = 0.75,
         cex.names=0.75)

```

Then, we can check if our data should be modeled using the Poisson distribution or Negative Binomial distribution. It can be helpful to plot the mean versus the variance of your data. Remember for the Poisson model, mean = variance, but for NB, mean < variance.

```{r include = TRUE, message=FALSE, warning=FALSE}
mean_counts <- apply(mRC_uni, 1, mean)
variance_counts <- apply(mRC_uni, 1, var)
df_cases <- data.frame(mean_counts, variance_counts)

ggplot(df_cases) +
        geom_point(aes(x=mean_counts, y=variance_counts)) + 
        geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
        scale_y_log10() +
        scale_x_log10()

```
Our data seems to fit a Negative Binomial distribution the most.


Running DESeqDataSetFromMatrix to generate the object class from deseq that stores the read counts and the intermediate estimated quantities during statistical analysis (`dds_uni`)

```{r include = TRUE, message=FALSE, warning=FALSE}

dds_uni <- DESeqDataSetFromMatrix(countData = mRC_uni,
                                  colData = sampleinfo_todosgroup,
                                  design = ~ batch + sex + age + group)

#Selecionar apenas miRNAs com pelo menos 18 amostras com no minimo 3 reads 
keep_counts <- rowSums(counts(dds_uni) >= 3) >= 18
dds_uni <- dds_uni[keep_counts,]

```

# Principal Component Analysis

```{r include = TRUE, message=FALSE, warning=FALSE} 
# Input is a matrix of log transformed values
rld <- rlog(dds_uni, blind=T)

plotPCA(rld, intgroup="group")
plotPCA(rld, intgroup="sex")
plotPCA(rld, intgroup="batch")
plotPCA(rld, intgroup="wave")

#NOTE: The plotPCA() function will only return the values for PC1 and PC2. If you would like to explore the additional PCs in your data or if you would like to identify genes that contribute most to the PCs, you can use the prcomp() function. For example, to plot any of the PCs we could run the following code:
  
# Input is a matrix of log transformed values
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
#grupo
df <- cbind(sampleinfo_todosgroup, pca$x)
PC1toPC5_var<- pca$sdev*pca$sdev
PC1toPC5_per <- round(PC1toPC5_var/sum(PC1toPC5_var)*100,1)
pc1 <- ggplot(df) + geom_point(aes(x=PC1, y=PC2, col = group), size = 3) +xlab(paste("PC1 - ", PC1toPC5_per[1], "%", sep = "")) + ylab(paste("PC2 - ", PC1toPC5_per[2], "%", sep = ""))
pc2 <- ggplot(df) + geom_point(aes(x=PC2, y=PC3, col = group), size = 3)+ xlab(paste("PC2 - ", PC1toPC5_per[2], "%", sep = "")) + ylab(paste("PC3 - ", PC1toPC5_per[3], "%", sep = ""))
library(ggpubr)
pc3 <- ggplot(df) + geom_point(aes(x=PC3, y=PC4, col = group), size = 3) +xlab(paste("PC3 - ", PC1toPC5_per[3], "%", sep = "")) + ylab(paste("PC4 - ", PC1toPC5_per[4], "%", sep = ""))
pc4 <- ggplot(df) + geom_point(aes(x=PC4, y=PC5, col = group), size = 3) +xlab(paste("PC4 - ", PC1toPC5_per[4], "%", sep = "")) + ylab(paste("PC5 - ", PC1toPC5_per[5], "%", sep = ""))

ggarrange(pc1, pc2, pc3, pc4,common.legend = TRUE,
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)

```

# Hierarchical Clustering

```{r include = TRUE, message=FALSE, warning=FALSE} 

library("pheatmap")

### Extract the rlog matrix from the object
rld_mat <- assay(rld)  

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

### Plot heatmap
df <- as.data.frame(colData(dds_uni)[,c("group","sex","wave")])
pheatmap(rld_cor, annotation_col= df)


##Heatmap##
select_uni <- order(rowMeans(counts(dds_uni,normalized=TRUE)),
                    decreasing=TRUE)[1:50]

df_uni <- as.data.frame(colData(dds_uni)[,c("sex","group")])

rownames(df_uni) <- gsub("[:MIMAT].*", "", rownames(df_uni))
rownames(ntd_uni) <- gsub("[:MIMAT].*", "", rownames(ntd_uni))
colnames(ntd_uni) <- gsub("[:Total].*","", colnames(ntd_uni))

pheatmap(assay(ntd_uni)[select_uni,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df_uni, fontsize_row=8, fontsize_col
         =8, main="Heatmap - top 50 expressed miRNA")

```
## DESeq2 differential gene expression analysis workflow
#Step 1: Estimate size factors

We can check size factors distribution

```{r include = TRUE, message=FALSE, warning=FALSE} 
sf <- sizeFactors(dds_uni)
ggplot() + aes(sf) + geom_histogram()
```


We can generate total normalized and non-normalized counts per sample and compare to size factors

```{r include = TRUE, message=FALSE, warning=FALSE} 
## Total number of non-normalized counts per sample
counts_per_sample <- colSums(counts(dds_uni))

## Total number of normalized counts per sample
norm_counts_per_sample <- colSums(counts(dds_uni, normalized=T))

ggplot(df) +
geom_point(aes(x=counts_per_sample, y=norm_counts_per_sample)) +
geom_line(aes(x=counts_per_sample, y=counts_per_sample, color="red"))

```

We can retrieve the normalized counts matrix from dds, we use the counts() function and add the argument normalized=TRUE.

```{r include = TRUE, message=FALSE, warning=FALSE} 
normalized_counts <- counts(dds_uni, normalized=TRUE)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

```
