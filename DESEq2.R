
csvfile <- file.path("trimmed_star", "Sample_tableAll.csv")
sampleTable <- read.csv(csvfile, row.names = 1)

sampleTable

filenames <- file.path("trimmed_star", paste0(sampleTable$Run, ".bam"))

file.exists(filenames)


library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)

seqinfo(bamfiles[1])

library("GenomicFeatures")

gtffile <- file.path("trimmed_star","Dioscorea_dumetorum_v1.0.gff")

(txdb <- makeTxDbFromGFF(gtffile, format="gff3", circ_seqs=character()))

(ebg <- exonsBy(txdb, by="gene"))


## Read counting step

library("GenomicAlignments")
library("BiocParallel")
library("tximeta")

register(SerialParam())

se <- summarizeOverlaps(features=ebg, reads=bamfiles,mode="Union",singleEnd=FALSE,ignore.strand=TRUE,fragments=TRUE )

dim(se)
head(row.names(se))

assayNames(se)

head(assay(se), 3)

colSums(assay(se))

rowRanges(se)

str(metadata(rowRanges(se)))

colData(se)

(colData(se) <- DataFrame(sampleTable))

se$SampleName

se$Conditions

## Test a different way

library("magrittr")

########################################################

## Exploratory analysis

se$Conditions <- relevel(se$Conditions, "4MAE")


round( colSums(assay(se)) / 1e6, 1 )

library("DESeq2")

dds <- DESeqDataSet(se, design = ~ SampleName + Conditions)

## First part to explore the data using normalize count

nrow(dds)
colData(dds)

## Pre-filtering

# Removing rows of the DESeqDataSet that have no counts, or only a single count across all samples.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
nrow(dds)


## (Effects of transformations on the variance)

#shifted logarithm transformation
# this gives log2(n + 1)

ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

##   The variance stabilizing transformation

vsd <- vst(dds, blind = FALSE)
meanSdPlot(assay(vsd))

## The rlog transformation

rld <- rlog(dds, blind = FALSE)
meanSdPlot(assay(rld))

library("dplyr")
library("ggplot2")


## Heatmap of the count matrix

dds <- estimateSizeFactors(dds)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Conditions","SampleName")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

## Heatmap with The variance stabilizing transformation

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

## Heatmap of the count matrix

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

## Heatmap of the sample-to-sample distances

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Conditions, vsd$SampleName, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## Principal component plot of the samples

plotPCA(vsd, intgroup=c("Conditions", "SampleName"))

#  PCA plot using the ggplot function

library("dplyr")
library("ggplot2")


pcaData <- plotPCA(vsd, intgroup=c("Conditions", "SampleName"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Conditions, shape=SampleName)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()



##########################################################################################################################

## Differential expression analysis

se$Conditions <- relevel(se$Conditions, "AH")

dds <- DESeqDataSet(se, design = ~ SampleName + Conditions)


## Pre-filtering

# A minimal pre-filtering to keep only rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
nrow(dds)

ddsAll <- DESeq(dds)

## 1. results 4DAH vs AH which is standard

resAll <- results(ddsAll)
head(resAll)
resAll

# the order of the coefficient as it appears in resultsNames(dds)

resultsNames(ddsAll)

#Log fold change shrinkage for visualization and ranking

resLFC4DAH_vs_AH <- lfcShrink(ddsAll, coef="Conditions_4DAH_vs_AH", type="apeglm")
resLFC4DAH_vs_AH

# we test for genes that show significant effects of treatment on gene counts more than doubling or less than halving, because 22=4

res4DAH_vs_AH <- results(ddsAll, lfcThreshold=2, alpha = 0.05)
res4DAH_vs_AH
summary(res4DAH_vs_AH)

#####################################################################

## 2. Results 3DAH vs AH

# the order of the coefficient as it appears in resultsNames(dds)

resultsNames(ddsAll)

#Log fold change shrinkage for visualization and ranking

resLFC3DAH_vs_AH <- lfcShrink(ddsAll, coef="Conditions_3DAH_vs_AH", type="apeglm")
resLFC3DAH_vs_AH

# we test for genes that show significant effects of treatment on gene counts more than doublingx2 or less than halving, because 22=4
# Used contrast to select 3DAH vs AH

res3DAH_vs_AH <- results(ddsAll, name="Conditions_3DAH_vs_AH", lfcThreshold=2, alpha = 0.05)
res3DAH_vs_AH 
summary(res3DAH_vs_AH)

#########################################################

## 3. Results 4MAE vs AH

# the order of the coefficient as it appears in resultsNames(dds)

resultsNames(ddsAll)

#Log fold change shrinkage for visualization and ranking

resLFC4MAE_vs_AH <- lfcShrink(ddsAll, coef="Conditions_4MAE_vs_AH", type="apeglm")
resLFC4MAE_vs_AH

# we test for genes that show significant effects of treatment 4MAE vs AH

res4MAE_vs_AH <- results(ddsAll, name="Conditions_4MAE_vs_AH", lfcThreshold=2, alpha = 0.05)
res4MAE_vs_AH 
summary(res4MAE_vs_AH)


##########################################################

## 4. Results 4DAH vs 3DAH

se$Conditions <- relevel(se$Conditions, "3DAH")

dds2 <- DESeqDataSet(se, design = ~ SampleName + Conditions)

## Pre-filtering

# A minimal pre-filtering to keep only rows that have at least 10 reads total
keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]
nrow(dds2)

ddsAll2 <- DESeq(dds2)

resAll2 <- results(ddsAll2)
head(resAll2)
resAll2

# the order of the coefficient as it appears in resultsNames(dds)

resultsNames(ddsAll2)

#Log fold change shrinkage for visualization and ranking

resLFC4DAH_vs_3DAH <- lfcShrink(ddsAll2, coef="Conditions_4DAH_vs_3DAH", type="apeglm")
resLFC4DAH_vs_3DAH

# 

res4DAH_vs_3DAH <- results(ddsAll2, lfcThreshold=2, alpha = 0.05)



#############################################################################################################################

## A. Genotype effect

#combine the factors of interest into a single factor with all combinations of the original factors

se$Conditions <- relevel(se$Conditions, "AH")

se$SampleName <- relevel(se$SampleName, "Ibosweet3")

ddsGeno_Cond1 <- DESeqDataSet(se, design = ~ SampleName + Conditions + SampleName:Conditions)


## Pre-filtering

# A minimal pre-filtering to keep only rows that have at least 10 reads total
keep <- rowSums(counts(ddsGeno_Cond1)) >= 10
ddsGeno_Cond1 <- ddsGeno_Cond1[keep,]
nrow(ddsGeno_Cond1)

ddsGeno_Cond1 <- DESeq(ddsGeno_Cond1)

resGeno_Cond1 <- results(ddsGeno_Cond1)
resGeno_Cond1
resultsNames(ddsGeno_Cond1)

##########################################

# A.1 the condition 4DAH_vs_AH effect for genotype Ibosweet3 (the main effect)

res_ibo_4DAH_vs_AH <- results(ddsGeno_Cond1, contrast=c("Conditions","4DAH","AH"), lfcThreshold=2, alpha = 0.05)

res_ibo_4DAH_vs_AH

###########################################

# A.1.1 the condition 3DAH_vs_AH effect for genotype Ibosweet3 (the main effect)


res_ibo_3DAH_vs_AH <- results(ddsGeno_Cond1, contrast=c("Conditions","3DAH","AH"), lfcThreshold=2, alpha = 0.05)

#############################################################################################

# A.1.1 the condition 4MAE_vs_AH effect for genotype Ibosweet3 (the main effect)


res_ibo_4MAE_vs_AH <- results(ddsGeno_Cond1, contrast=c("Conditions","4MAE","AH"), lfcThreshold=2, alpha = 0.05)


###########################################################################################

## A.2 the condition 4DAH_vs_AH effect for genotype Bayangam2 (the main effect)
# the condition effect for Bayangam2.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Bayangam2 compared to Ibosweet3)

res_baya_4DAH_vs_AH <- results(ddsGeno_Cond1, contrast=list( c("Conditions_4DAH_vs_AH","SampleNameBayangam2.Conditions4DAH") ), lfcThreshold=2, alpha = 0.05)

########################################################################

## A.2.2 the condition 3DAH_vs_AH effect for genotype Bayangam2 (the main effect)
# the condition effect for Bayangam2.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Bayangam2 compared to Ibosweet3)

res_baya_3DAH_vs_AH <- results(ddsGeno_Cond1, contrast=list( c("Conditions_3DAH_vs_AH","SampleNameBayangam2.Conditions3DAH") ), lfcThreshold=2, alpha = 0.05)

################################################################################

## A.2.2 the condition 4MAE_vs_AH effect for genotype Bayangam2 (the main effect)
# the condition effect for Bayangam2.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Bayangam2 compared to Ibosweet3)

res_baya_4MAE_vs_AH <- results(ddsGeno_Cond1, contrast=list( c("Conditions_4MAE_vs_AH","SampleNameBayangam2.Conditions4MAE") ), lfcThreshold=2, alpha = 0.05)

############################################################################

## A.3 the condition 4DAH_vs_AH effect for genotype Bangou1 (the main effect)
# the condition effect for Bangou1.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Bangou1 compared to Ibosweet3)

res_bang_4DAH_vs_AH <- results(ddsGeno_Cond1, contrast=list( c("Conditions_4DAH_vs_AH","SampleNameBangou1.Conditions4DAH") ), lfcThreshold=2, alpha = 0.05)

############################################################################


## A.3.2 the condition 3DAH_vs_AH effect for genotype Bangou1 (the main effect)
# the condition effect for Bangou1.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Bangou1 compared to Ibosweet3)

res_bang_3DAH_vs_AH <- results(ddsGeno_Cond1, contrast=list( c("Conditions_3DAH_vs_AH","SampleNameBangou1.Conditions3DAH") ), lfcThreshold=2, alpha = 0.05)

#################################################################################

## A.3.2 the condition 4MAE_vs_AH effect for genotype Bangou1 (the main effect)
# the condition effect for Bangou1.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Bangou1 compared to Ibosweet3)

res_bang_4MAE_vs_AH <- results(ddsGeno_Cond1, contrast=list( c("Conditions_4MAE_vs_AH","SampleNameBangou1.Conditions4MAE") ), lfcThreshold=2, alpha = 0.05)

################################################################################

## A.4 the condition 4DAH_vs_AH effect for genotype Fonkouankem1 (the main effect)
# the condition effect for Fonkouankem1.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Fonkouankem1 compared to Ibosweet3)

res_fonk_4DAH_vs_AH <- results(ddsGeno_Cond1, contrast=list( c("Conditions_4DAH_vs_AH","SampleNameFonkouankem1.Conditions4DAH") ), lfcThreshold=2, alpha = 0.05)


####################################################################################

## A.4 the condition 3DAH_vs_AH effect for genotype Fonkouankem1 (the main effect)
# the condition effect for Fonkouankem1.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Fonkouankem1 compared to Ibosweet3)

res_fonk_3DAH_vs_AH <- results(ddsGeno_Cond1, contrast=list( c("Conditions_3DAH_vs_AH","SampleNameFonkouankem1.Conditions3DAH") ), lfcThreshold=2, alpha = 0.05)

#############################################################################

## A.4 the condition 4MAE_vs_AH effect for genotype Fonkouankem1 (the main effect)
# the condition effect for Fonkouankem1.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Fonkouankem1 compared to Ibosweet3)

res_fonk_4MAE_vs_AH <- results(ddsGeno_Cond1, contrast=list( c("Conditions_4MAE_vs_AH","SampleNameFonkouankem1.Conditions4MAE") ), lfcThreshold=2, alpha = 0.05)


#############################################################################

## B. Genotype effect plus 3DAH as base of comparison


#combine the factors of interest into a single factor with all combinations of the original factors

se$Conditions <- relevel(se$Conditions, "3DAH")

se$SampleName <- relevel(se$SampleName, "Ibosweet3")

ddsGeno_Cond3DAH <- DESeqDataSet(se, design = ~ SampleName + Conditions + SampleName:Conditions)


## Pre-filtering

# A minimal pre-filtering to keep only rows that have at least 10 reads total
keep <- rowSums(counts(ddsGeno_Cond3DAH)) >= 10
ddsGeno_Cond3DAH <- ddsGeno_Cond3DAH[keep,]
nrow(ddsGeno_Cond3DAH)

ddsGeno_Cond3DAH <- DESeq(ddsGeno_Cond3DAH)

resGeno_Cond3DAH <- results(ddsGeno_Cond3DAH)
resGeno_Cond3DAH

resultsNames(ddsGeno_Cond3DAH)

##############################################


# B.1 the condition 4DAH_vs_3DAH effect for genotype Ibosweet3 (the main effect)

res_ibo_4DAH_vs_3DAH <- results(ddsGeno_Cond3DAH, contrast=c("Conditions","4DAH","3DAH"), lfcThreshold=2, alpha = 0.05)

#####################################################################


## B.2 the condition 4DAH_vs_3DAH effect for genotype Bayangam2 (the main effect)
# the condition effect for Bayangam2.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Bayangam2 compared to Ibosweet3)

res_baya_4DAH_vs_3DAH <- results(ddsGeno_Cond3DAH, contrast=list( c("Conditions_4DAH_vs_3DAH","SampleNameBayangam2.Conditions4DAH") ), lfcThreshold=2, alpha = 0.05)


###########################################################################


## B.3 the condition 4DAH_vs_3DAH effect for genotype Bangou1 (the main effect)
# the condition effect for Bangou1.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Bangou1 compared to Ibosweet3)

res_bang_4DAH_vs_3DAH <- results(ddsGeno_Cond3DAH, contrast=list( c("Conditions_4DAH_vs_3DAH","SampleNameBangou1.Conditions4DAH") ), lfcThreshold=2, alpha = 0.05)

############################################################################

## B.4 the condition 4DAH_vs_3DAH effect for genotype Fonkouankem1 (the main effect)
# the condition effect for Fonkouankem1.
# this is the main effect *plus* the interaction term
# (the extra condition effect in Fonkouankem1 compared to Ibosweet3)

res_fonk_4DAH_vs_3DAH <- results(ddsGeno_Cond3DAH, contrast=list( c("Conditions_4DAH_vs_3DAH","SampleNameFonkouankem1.Conditions4DAH") ), lfcThreshold=2, alpha = 0.05)


###################################################################################################################################

## Time course experiments

library("DESeq2")

se$Conditions <- relevel(se$Conditions, "AH")

se$SampleName <- relevel(se$SampleName, "Ibosweet3")


ddsTC <- DESeqDataSet(se, ~ SampleName + Conditions + SampleName:Conditions)



keep <- rowSums(counts(ddsTC)) >= 10
ddsTC <- ddsTC[keep,]
nrow(ddsTC)


ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ SampleName + Conditions)

resTC <- results(ddsTC)

resTC$symbol <- mcols(ddsTC)$symbol

head(resTC[order(resTC$padj),], 4)

fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("Conditions","SampleName"), returnData = TRUE)

library(ggplot2)
fiss$Conditions <- as.factor(fiss$Conditions)
ggplot(fiss,
       aes(x = Conditions, y = count, color = SampleName, group = SampleName)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10()

resultsNames(ddsTC)

## SampleNameBayangam2.Conditions3DAH

resBaya3DAHTC <- results(ddsTC, name="SampleNameBayangam2.Conditions3DAH", test="Wald", lfcThreshold=2, alpha = 0.05)

#SampleNameBayangam2.Conditions4DAH

resBaya4DAHTC <- results(ddsTC, name="SampleNameBayangam2.Conditions4DAH", test="Wald", lfcThreshold=2, alpha = 0.05)


#SampleNameBayangam2.Conditions4MAE

resBaya4MAETC <- results(ddsTC, name="SampleNameBayangam2.Conditions4MAE", test="Wald", lfcThreshold=2, alpha = 0.05)

#SampleNameBangou1.Conditions3DAH

resBang3DAHTC <- results(ddsTC, name="SampleNameBangou1.Conditions3DAH", test="Wald", lfcThreshold=2, alpha = 0.05)

#SampleNameBangou1.Conditions4DAH

resBang4DAHTC <- results(ddsTC, name="SampleNameBangou1.Conditions4DAH", test="Wald", lfcThreshold=2, alpha = 0.05)


#SampleNameBangou1.Conditions4MAE

resBang4MAETC <- results(ddsTC, name="SampleNameBangou1.Conditions4MAE", test="Wald", lfcThreshold=2, alpha = 0.05)

#SampleNameFonkouankem1.Conditions3DAH

resFonk3DAHTC <- results(ddsTC, name="SampleNameFonkouankem1.Conditions3DAH", test="Wald", lfcThreshold=2, alpha = 0.05)

#SampleNameFonkouankem1.Conditions4DAH

resFonk4DAHTC <- results(ddsTC, name="SampleNameFonkouankem1.Conditions4DAH", test="Wald", lfcThreshold=2, alpha = 0.05)

#SampleNameFonkouankem1.Conditions4MAE

resFonk4MAETC <- results(ddsTC, name="SampleNameFonkouankem1.Conditions4MAE", test="Wald", lfcThreshold=2, alpha = 0.05)


# Ibo Conditions_4MAE_vs_AH

res_ibo__4MAE_vs_AHTC <- results(ddsTC, name="Conditions_4MAE_vs_AH", test="Wald", lfcThreshold=2, alpha = 0.05)


## Ibo Conditions_3DAH_vs_AH

res_ibo__3DAH_vs_AHTC <- results(ddsTC, name="Conditions_3DAH_vs_AH", test="Wald", lfcThreshold=2, alpha = 0.05)


# Ibo Conditions_4DAH_vs_AH

res_ibo__4DAH_vs_AHTC <- results(ddsTC, name="Conditions_4DAH_vs_AH", test="Wald", lfcThreshold=2, alpha = 0.05)

#SampleName_Bangou1_vs_Ibosweet3

res_Bang_vs_IboTC <- results(ddsTC, name="SampleName_Bangou1_vs_Ibosweet3", test="Wald", lfcThreshold=2, alpha = 0.05)

#SampleName_Fonkouankem1_vs_Ibosweet3

res_Fonk_vs_IboTC <- results(ddsTC, name="SampleName_Fonkouankem1_vs_Ibosweet3", test="Wald", lfcThreshold=2, alpha = 0.05)

## SampleName_Bayangam2_vs_Ibosweet3

res_Baya_vs_IboTC <- results(ddsTC, name="SampleName_Bayangam2_vs_Ibosweet3", test="Wald", lfcThreshold=2, alpha = 0.05)

