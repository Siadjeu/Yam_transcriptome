## Load data

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


##Read counting step

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

##(Reference: https://f1000research.com/articles/5-1438/v2)
## effect conditions no matter the sample name

group <- paste(colData(se)$Conditions, sep=".")
group <- factor(group)

group 
table(group)
levels(group)

## Change factor name because The levels must by syntactically valid names in R, see help(make.names)
levels(group)[c(1,2,3)] <- c('DAH3','DAH4','MAE4')
group
table(group)
levels(group)

# DGEList
library(edgeR)

edger <- DGEList(assays(se)$counts, group=group)

edger$samples

names(edger)

# Filtering to remove low counts
# For the current analysis, we keep genes that have at least 10 counts

keep <- rowSums(cpm(edger) > 0.5) >= 2
table(keep)

## because the DGEList object is subsetted to retain only the non-filtered genes we run this code
y <- edger[keep, , keep.lib.sizes=FALSE]
y


y <- calcNormFactors(y)

y$samples

pch <- c(0,1,2,5,6,7,15,16,17,18,19,21,22,23,24,25)
colors <- rep(c("green", "red", "blue", "black"), 12)
plotMDS(y, col=colors[group], pch=pch[group])
legend("bottomleft", legend=levels(group), pch=pch, col=colors, ncol=2)

plotMD(y, column=11)
abline(h=0, col="red", lty=2, lwd=2)

##Design matrix

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

##Dispersion estimation  
library(statmod)
y <- estimateDisp(y, design, robust=TRUE)

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

plotQLDisp(fit)

summary(fit$df.prior)

##Differential expression analysis

# comparisons is that between the AH and DAH3

trDAH3vsAH <- makeContrasts(DAH3-AH, levels=design)

tr1DAH3vsAH <- glmTreat(fit, contrast=trDAH3vsAH, lfc=2)

EdgeDAH3vsAH.52 <- topTags(tr1DAH3vsAH, n = 52)
EdgeDAH3vsAH.300 <- topTags(tr1DAH3vsAH, n = 300)


##############

# comparisons is that between the AH and DAH4

trDAH4vsAH <- makeContrasts(DAH4-AH, levels=design)

tr1DAH4vsAH <- glmTreat(fit, contrast=trDAH4vsAH, lfc=2)

EdgeDAH4vsAH.11 <- topTags(tr1DAH4vsAH, n = 11)
EdgeDAH4vsAH.300 <- topTags(tr1DAH4vsAH, n = 300)


is.deDAH4vsAH <- decideTestsDGE(tr1DAH4vsAH)
summary(is.deDAH4vsAH)


#######
# comparisons is that between the AH and MAE4

trMAE4vsAH <- makeContrasts(MAE4-AH, levels=design)

tr1MAE4vsAH <- glmTreat(fit, contrast=trMAE4vsAH, lfc=2)

EdgeMAE4vsAH.1617 <- topTags(tr1MAE4vsAH, n = 1617)
EdgeMAE4vsAH.2000 <- topTags(tr1MAE4vsAH, n = 2000)


is.deMAE4vsAH <- decideTestsDGE(tr1MAE4vsAH)
summary(is.deMAE4vsAH)



##########################################################################################################################

group <- paste(colData(se)$Conditions, colData(se)$SampleName, sep=".")

group <- factor(group)

group 
table(group)
## change factor name because The levels must by syntactically valid names in R, see help(make.names)
levels(group)[c(1,2,3,4,5,6,7,8,9,10,11,12)] <- c('DAH3.Bangou1','DAH3.Bayangam2','DAH3.Fonkouankem1','DAH3.Ibosweet3',
                                                  'DAH4.Bangou1','DAH4.Bayangam2','DAH4.Fonkouankem1','DAH4.Ibosweet3','MAE4.Bangou1',
                                                  'MAE4.Bayangam2','MAE4.Fonkouankem1','MAE4.Ibosweet3')
group
table(group)
levels(group)
edger <- DGEList(assays(se)$counts, group=group)



edger$samples

names(edger)

# Filtering to remove low counts
# For the current analysis, we keep genes that have at least 10 counts

##keep <-rowSums(edger$counts) >=10 it is not the truth

##table(keep)

## because the DGEList object is subsetted to retain only the non-filtered genes we run this code
y <- edger[keep, , keep.lib.sizes=FALSE]
y



##Normalization for composition bias
y <- calcNormFactors(y)

y$samples

pch <- c(0,1,2,5,6,7,15,16,17,18,19,21,22,23,24,25)
colors <- rep(c("darkgreen", "red", "blue", "black"), 3)
plotMDS(y, col=colors[group], pch=pch[group])
legend("bottomleft", legend=levels(group), pch=pch, col=colors, ncol=2)

plotMD(y, column=11)
abline(h=0, col="red", lty=2, lwd=2)

##Design matrix

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

##Dispersion estimation  
library(statmod)
y <- estimateDisp(y, design, robust=TRUE)

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

plotQLDisp(fit)

summary(fit$df.prior)

#####################################################################################################################


##Differential expression analysis

# Bayangam: comparisons is that between the AH.Bayangam2 and DAH3.Bayangam2

trDAH3.BayavsAH.Baya <- makeContrasts(DAH3.Bayangam2-AH.Bayangam2, levels=design)

tr1DAH3.BayavsAH.Baya <- glmTreat(fit, contrast=trDAH3.BayavsAH.Baya, lfc=2)

EdgeDAH3.BayavsAH.Baya20 <- topTags(tr1DAH3.BayavsAH.Baya, n = 20)
EdgeDAH3.BayavsAH.Baya300 <- topTags(tr1DAH3.BayavsAH.Baya, n = 300)


## Bayangam: comparisons is that between the AH.Bayangam2 and DAH4.Bayangam2


trDAH4.BayavsAH.Baya <- makeContrasts(DAH4.Bayangam2-AH.Bayangam2, levels=design)

tr1DAH4.BayavsAH.Baya <- glmTreat(fit, contrast=trDAH4.BayavsAH.Baya, lfc=2)

EdgeDAH4.BayavsAH.Baya20 <- topTags(tr1DAH4.BayavsAH.Baya, n = 20)
EdgeDAH4.BayavsAH.Baya300 <- topTags(tr1DAH4.BayavsAH.Baya, n = 300)


## Bayangam: comparisons is that between the DAH3.Bayangam2 and DAH4.Bayangam2

trDAH4.BayavsDAH3.Baya <- makeContrasts(DAH4.Bayangam2-DAH3.Bayangam2, levels=design)

tr1DAH4.BayavsDAH3.Baya <- glmTreat(fit, contrast=trDAH4.BayavsDAH3.Baya, lfc=2)

EdgeDAH4.BayavsDAH3.Baya20 <- topTags(tr1DAH4.BayavsDAH3.Baya, n = 20)
EdgeDAH4.BayavsDAH3.Baya300 <- topTags(tr1DAH4.BayavsDAH3.Baya, n = 300)


## Bayangam: comparisons is that between the AH.Bayangam2 and MAE4.Bayangam2

trMAE4.BayavsAH.Baya <- makeContrasts(MAE4.Bayangam2-AH.Bayangam2, levels=design)

tr1MAE4.BayavsAH.Baya <- glmTreat(fit, contrast=trMAE4.BayavsAH.Baya, lfc=2)

EdgeMAE4.BayavsAH.Baya20 <- topTags(tr1MAE4.BayavsAH.Baya, n = 20)
EdgeMAE4.BayavsAH.Baya300 <- topTags(tr1MAE4.BayavsAH.Baya, n = 300)



## Bayangam: comparisons is that between the MAE4.Bayangam2 and AH.Bayangam2

trAH.BayavsMAE4.Baya <- makeContrasts(AH.Bayangam2-MAE4.Bayangam2, levels=design)

tr1AH.BayavsMAE4.Baya <- glmTreat(fit, contrast=trAH.BayavsMAE4.Baya, lfc=2)

EdgeAH.BayavsMAE4.Baya20 <- topTags(tr1AH.BayavsMAE4.Baya, n = 20)
EdgeAH.BayavsMAE4.Baya300 <- topTags(tr1AH.BayavsMAE4.Baya, n = 300)



#######################################################################


# Fonkouankem1: comparisons is that between the AH.Fonkouankem1 and DAH3.Fonkouankem1

trDAH3.FonkvsAH.Fonk <- makeContrasts(DAH3.Fonkouankem1-AH.Fonkouankem1, levels=design)

tr1DAH3.FonkvsAH.Fonk <- glmTreat(fit, contrast=trDAH3.FonkvsAH.Fonk, lfc=2)

EdgeDAH3.FonkvsAH.Fonk20 <- topTags(tr1DAH3.FonkvsAH.Fonk, n = 20)
EdgeDAH3.FonkvsAH.Fonk300 <- topTags(tr1DAH3.FonkvsAH.Fonk, n = 300)


is.deDAH3.FonkvsAH.Fonk <- decideTestsDGE(tr1DAH3.FonkvsAH.Fonk)
summary(is.deDAH3.FonkvsAH.Fonk)


## Fonkouankem1: comparisons is that between the AH.Fonkouankem1 and DAH4.Fonkouankem1


trDAH4.FonkvsAH.Fonk <- makeContrasts(DAH4.Fonkouankem1-AH.Fonkouankem1, levels=design)

tr1DAH4.FonkvsAH.Fonk <- glmTreat(fit, contrast=trDAH4.FonkvsAH.Fonk, lfc=2)

EdgeDAH4.FonkvsAH.Fonk20 <- topTags(tr1DAH4.FonkvsAH.Fonk, n = 20)
EdgeDAH4.FonkvsAH.Fonk300 <- topTags(tr1DAH4.FonkvsAH.Fonk, n = 300)

## Fonkouankem1: comparisons is that between the DAH3.Fonkouankem1 and DAH4.Fonkouankem1

trDAH4.FonkvsDAH3.Fonk <- makeContrasts(DAH4.Fonkouankem1-DAH3.Fonkouankem1, levels=design)

tr1DAH4.FonkvsDAH3.Fonk <- glmTreat(fit, contrast=trDAH4.FonkvsDAH3.Fonk, lfc=2)

EdgeDAH4.FonkvsDAH3.Fonk20 <- topTags(tr1DAH4.FonkvsDAH3.Fonk, n = 20)
EdgeDAH4.FonkvsDAH3.Fonk300 <- topTags(tr1DAH4.FonkvsDAH3.Fonk, n = 300)


## Fonkngam: comparisons is that between the MAE4.Fonkouankem1 and AH.Fonkouankem1

trAH.FonkvsMAE4.Fonk <- makeContrasts(AH.Fonkouankem1-MAE4.Fonkouankem1, levels=design)

tr1AH.FonkvsMAE4.Fonk <- glmTreat(fit, contrast=trAH.FonkvsMAE4.Fonk, lfc=2)

EdgeAH.FonkvsMAE4.Fonk20 <- topTags(tr1AH.FonkvsMAE4.Fonk, n = 20)
EdgeAH.FonkvsMAE4.Fonk300 <- topTags(tr1AH.FonkvsMAE4.Fonk, n = 300)


###################################################

# Bangou1: comparisons is that between the AH.Bangou1 and DAH3.Bangou1

trDAH3.BangvsAH.Bang <- makeContrasts(DAH3.Bangou1-AH.Bangou1, levels=design)

tr1DAH3.BangvsAH.Bang <- glmTreat(fit, contrast=trDAH3.BangvsAH.Bang, lfc=2)

EdgeDAH3.BangvsAH.Bang20 <- topTags(tr1DAH3.BangvsAH.Bang, n = 20)
EdgeDAH3.BangvsAH.Bang300 <- topTags(tr1DAH3.BangvsAH.Bang, n = 300)


## Bangou1: comparisons is that between the AH.Bangou1 and DAH4.Bangou1


trDAH4.BangvsAH.Bang <- makeContrasts(DAH4.Bangou1-AH.Bangou1, levels=design)

tr1DAH4.BangvsAH.Bang <- glmTreat(fit, contrast=trDAH4.BangvsAH.Bang, lfc=2)

EdgeDAH4.BangvsAH.Bang20 <- topTags(tr1DAH4.BangvsAH.Bang, n = 20)
EdgeDAH4.BangvsAH.Bang300 <- topTags(tr1DAH4.BangvsAH.Bang, n = 300)


## Bangou1: comparisons is that between the DAH3.Bangou1 and DAH4.Bangou1

trDAH4.BangvsDAH3.Bang <- makeContrasts(DAH4.Bangou1-DAH3.Bangou1, levels=design)

tr1DAH4.BangvsDAH3.Bang <- glmTreat(fit, contrast=trDAH4.BangvsDAH3.Bang, lfc=2)

EdgeDAH4.BangvsDAH3.Bang20 <- topTags(tr1DAH4.BangvsDAH3.Bang, n = 20)
EdgeDAH4.BangvsDAH3.Bang300 <- topTags(tr1DAH4.BangvsDAH3.Bang, n = 300)


is.deDAH4.BangvsDAH3.Bang <- decideTestsDGE(tr1DAH4.BangvsDAH3.Bang)
summary(is.deDAH4.BangvsDAH3.Bang)

# export the subset top 20 and top 300

write.csv(as.data.frame(EdgeAH.BangvsMAE4.Bang20), 
          file="Edge_AH.BangvsMAE4.Bangtop20[2-0.05].csv")


write.csv(as.data.frame(EdgeAH.BangvsMAE4.Bang300), 
          file="Edge_AH.BangvsMAE4.Bangtop300[2-0.05].csv")

#######################################################################################################
# Ibosweet3: comparisons is that between the AH.Ibosweet3 and DAH3.Ibosweet3

trDAH3.IbovsAH.Ibo <- makeContrasts(DAH3.Ibosweet3-AH.Ibosweet3, levels=design)

tr1DAH3.IbovsAH.Ibo <- glmTreat(fit, contrast=trDAH3.IbovsAH.Ibo, lfc=2)

EdgeDAH3.IbovsAH.Ibo20 <- topTags(tr1DAH3.IbovsAH.Ibo, n = 20)
EdgeDAH3.IbovsAH.Ibo300 <- topTags(tr1DAH3.IbovsAH.Ibo, n = 300)


## Ibosweet3: comparisons is that between the AH.Ibosweet3 and DAH4.Ibosweet3


trDAH4.IbovsAH.Ibo <- makeContrasts(DAH4.Ibosweet3-AH.Ibosweet3, levels=design)

tr1DAH4.IbovsAH.Ibo <- glmTreat(fit, contrast=trDAH4.IbovsAH.Ibo, lfc=2)

EdgeDAH4.IbovsAH.Ibo20 <- topTags(tr1DAH4.IbovsAH.Ibo, n = 20)
EdgeDAH4.IbovsAH.Ibo300 <- topTags(tr1DAH4.IbovsAH.Ibo, n = 300)


## Ibosweet3: comparisons is that between the DAH3.Ibosweet3 and DAH4.Ibosweet3

trDAH4.IbovsDAH3.Ibo <- makeContrasts(DAH4.Ibosweet3-DAH3.Ibosweet3, levels=design)

tr1DAH4.IbovsDAH3.Ibo <- glmTreat(fit, contrast=trDAH4.IbovsDAH3.Ibo, lfc=2)

EdgeDAH4.IbovsDAH3.Ibo20 <- topTags(tr1DAH4.IbovsDAH3.Ibo, n = 20)
EdgeDAH4.IbovsDAH3.Ibo300 <- topTags(tr1DAH4.IbovsDAH3.Ibo, n = 300)



## Ibo: comparisons is that between the MAE4.Ibosweet3 and AH.Ibosweet3

trAH.IbovsMAE4.Ibo <- makeContrasts(AH.Ibosweet3-MAE4.Ibosweet3, levels=design)

tr1AH.IbovsMAE4.Ibo <- glmTreat(fit, contrast=trAH.IbovsMAE4.Ibo, lfc=2)

EdgeAH.IbovsMAE4.Ibo20 <- topTags(tr1AH.IbovsMAE4.Ibo, n = 20)
EdgeAH.IbovsMAE4.Ibo300 <- topTags(tr1AH.IbovsMAE4.Ibo, n = 300)


########################################################################################################

## Interaction Ibo and bayangam for AH and DAH3

IntBaya_Ibo3DAH_vs_AH <- makeContrasts((DAH3.Bayangam2-DAH3.Ibosweet3)-(AH.Bayangam2-AH.Ibosweet3), levels=design)

resIntBay_Ibo3DAH_vs_AH <- glmTreat(fit, contrast=IntBaya_Ibo3DAH_vs_AH, lfc=2)

Edge_interAH.IbovsDAH3.Baya20 <- topTags(resIntBay_Ibo3DAH_vs_AH, n = 20)
Edge_interAH.IbovsDAH3.Baya300 <- topTags(resIntBay_Ibo3DAH_vs_AH, n = 300)


## Interaction Ibo and Bayangam2 for AH and DAH4

IntBaya_Ibo4DAH_vs_AH <- makeContrasts((DAH4.Bayangam2-DAH4.Ibosweet3)-(AH.Bayangam2-AH.Ibosweet3), levels=design)

resIntBay_Ibo4DAH_vs_AH <- glmTreat(fit, contrast=IntBaya_Ibo4DAH_vs_AH, lfc=2)

Edge_interAH.IbovsDAH4.Baya20 <- topTags(resIntBay_Ibo4DAH_vs_AH, n = 20)
Edge_interAH.IbovsDAH4.Baya300 <- topTags(resIntBay_Ibo4DAH_vs_AH, n = 300)


##############################################

## Interaction Ibosweet3) and Fonkouankem1 for AH and DAH3

IntFonka_Ibo3DAH_vs_AH <- makeContrasts((DAH3.Fonkouankem1-DAH3.Ibosweet3)-(AH.Fonkouankem1-AH.Ibosweet3), levels=design)

resIntFonk_Ibo3DAH_vs_AH <- glmTreat(fit, contrast=IntFonka_Ibo3DAH_vs_AH, lfc=2)

Edge_interAH.IbovsDAH3.Fonka20 <- topTags(resIntFonk_Ibo3DAH_vs_AH, n = 20)
Edge_interAH.IbovsDAH3.Fonka300 <- topTags(resIntFonk_Ibo3DAH_vs_AH, n = 300)

is.de_interAH.IbovsDAH3.Fonka20 <- decideTestsDGE(resIntFonk_Ibo3DAH_vs_AH)
summary(is.de_interAH.IbovsDAH3.Fonka20)


## Interaction Ibosweet3) and Fonkouankem1 for AH and DAH4

IntFonka_Ibo4DAH_vs_AH <- makeContrasts((DAH4.Fonkouankem1-DAH4.Ibosweet3)-(AH.Fonkouankem1-AH.Ibosweet3), levels=design)

resIntFonk_Ibo4DAH_vs_AH <- glmTreat(fit, contrast=IntFonka_Ibo4DAH_vs_AH, lfc=2)

Edge_interAH.IbovsDAH4.Fonka20 <- topTags(resIntFonk_Ibo4DAH_vs_AH, n = 20)
Edge_interAH.IbovsDAH4.Fonka300 <- topTags(resIntFonk_Ibo4DAH_vs_AH, n = 300)


####################################

## Interaction Ibosweet3) and Bangou1 for AH and DAH3

IntBanga_Ibo3DAH_vs_AH <- makeContrasts((DAH3.Bangou1-DAH3.Ibosweet3)-(AH.Bangou1-AH.Ibosweet3), levels=design)

resIntBang_Ibo3DAH_vs_AH <- glmTreat(fit, contrast=IntBanga_Ibo3DAH_vs_AH, lfc=2)

Edge_interAH.IbovsDAH3.Banga20 <- topTags(resIntBang_Ibo3DAH_vs_AH, n = 20)
Edge_interAH.IbovsDAH3.Banga300 <- topTags(resIntBang_Ibo3DAH_vs_AH, n = 300)


## Interaction Ibosweet3) and Bangou1 for AH and DAH4

IntBanga_Ibo4DAH_vs_AH <- makeContrasts((DAH4.Bangou1-DAH4.Ibosweet3)-(AH.Bangou1-AH.Ibosweet3), levels=design)

resIntBang_Ibo4DAH_vs_AH <- glmTreat(fit, contrast=IntBanga_Ibo4DAH_vs_AH, lfc=2)

Edge_interAH.IbovsDAH4.Banga20 <- topTags(resIntBang_Ibo4DAH_vs_AH, n = 20)
Edge_interAH.IbovsDAH4.Banga300 <- topTags(resIntBang_Ibo4DAH_vs_AH, n = 300)


########################################################################################################

