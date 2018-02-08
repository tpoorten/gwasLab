################
## LAB - Computing Solutions: Genome-Wide Association Study (GWAS) 

## Instructor: Tom Poorten

setwd("~/gwasLab/")
# setwd("~/Downloads/gwasLab/")
rm(list = ls())

# source("https://bioconductor.org/biocLite.R")
# biocLite("SNPRelate")
# biocLite("GWASTools")

library(SNPRelate)
# https://bioconductor.org/packages/release/bioc/html/SNPRelate.html
library(GWASTools)
# https://bioconductor.org/packages/release/bioc/html/GWASTools.html
#
library(ggplot2)
# //////////////

################
## Prepare SNP dataset

# convert Plink file to GDS file
snpgdsBED2GDS(bed.fn = "genotype_data/HDRA-G6-4-RDP1-filter2.bed",
              bim.fn = "genotype_data/HDRA-G6-4-RDP1-filter2.bim", 
              fam.fn = "genotype_data/HDRA-G6-4-RDP1-filter2.fam", 
              out.gdsfn =  "genotype_data/HDRA-G6-4-RDP1-filter2.gds")

# open GDS file, then genofile is used to access it; be sure to close with snpgdsClose(genofile)
genofile <- snpgdsOpen("genotype_data/HDRA-G6-4-RDP1-filter2.gds", readonly = FALSE)

geno.sample.ids = read.gdsn(index.gdsn(genofile, "sample.id"))
head(geno.sample.ids)
length(geno.sample.ids)

## Quick QC
sampleMissingRate <- snpgdsSampMissRate(genofile)
summary(sampleMissingRate)
hist(sampleMissingRate, breaks = 30)
rm(sampleMissingRate)

snpRateFreqs <- snpgdsSNPRateFreq(genofile, with.snp.id=TRUE)
head(data.frame(snpRateFreqs))
summary(snpRateFreqs$MissingRate)
hist(snpRateFreqs$MissingRate)
hist(snpRateFreqs$AlleleFreq)
hist(snpRateFreqs$AlleleFreq, breaks=__)
rm(snpRateFreqs)
# //////////////

################
## Prepare Phenotype dataset
sampleInfo = read.table(file = "phenotype_data/phenoAvLen_G6_4_RDP12_ALL.txt", header=T, stringsAsFactors = F)
IDtab = read.table("phenotype_data/HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv", header=F, sep="\t", stringsAsFactors = F, quote="")

head(sampleInfo)
head(IDtab)
sampleInfo$ID = IDtab$__[match(sampleInfo$IID, IDtab$V3)]
sampleInfo$pop = IDtab$__[match(sampleInfo$IID, IDtab$V3)]
# subset and sort, based on geno.sample.ids
sampleInfo = sampleInfo[match(geno.sample.ids , sampleInfo$ID),]
sampleInfo = na.omit(sampleInfo)

hist(sampleInfo$AVERAGE_LENGTH)
hist(sampleInfo$AVERAGE_LENGTH, breaks=__)
sampleInfo$AVERAGE_LENGTH2 = ifelse(sampleInfo$AVERAGE_LENGTH > __, 1, 0)

# reset to ensure that list is updated 
geno.sample.ids = sampleInfo$ID
rm(IDtab)
# //////////////

################
## Population structure analysis

## LD pruning (SNPRelate) 
snpset.id0 <- snpgdsSelectSNP(genofile, missing.rate=0.05, maf = 0.05, sample.id = geno.sample.ids)
length(snpset.id0)
set.seed(1000) ; snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, slide.max.bp = 500000, sample.id = geno.sample.ids, snp.id = snpset.id0)

# Get selected snp id
snpset.id <- unlist(snpset)

## Run PCA
pcaLDP <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, sample.id = geno.sample.ids)

# variance proportion (%)
pc.percentLDP <- pcaLDP$varprop*100
head(round(pc.percentLDP, 2))
barplot(pc.percentLDP[1:15])

## PCA plots
identical(sampleInfo$ID, pcaLDP$sample.id) # check that sample ids match up
sampleInfo$PC1LDP = pcaLDP$eigenvect[,__]
sampleInfo$PC2LDP = pcaLDP$eigenvect[,__]
sampleInfo$PC3LDP = pcaLDP$eigenvect[,__]

ggplot(sampleInfo, aes(PC1LDP,PC2LDP))  + 
  theme_bw() +
  geom_point(alpha=.8, size=1, aes(color=pop)) + 
  xlab(paste0("PC1 (",round(pc.percentLDP[1],1),"%)")) +
  ylab(paste0("PC2 (",round(pc.percentLDP[2],1),"%)")) +
  coord_equal()
#
ggplot(sampleInfo, aes(PC1LDP,PC3LDP))  + 
  theme_bw() +
  geom_point(alpha=.8, size=1, aes(color=pop)) + 
  xlab(paste0("PC1 (",round(pc.percentLDP[1],1),"%)")) +
  ylab(paste0("PC3 (",round(pc.percentLDP[3],1),"%)")) +
  coord_equal()

# PC vs phenotype
ggplot(sampleInfo, aes(PC1LDP, AVERAGE_LENGTH, color=pop)) + geom_point()
ggplot(sampleInfo, aes(PC2LDP, AVERAGE_LENGTH, color=pop)) + geom_point()
# //////////////

################
## Relatedness analysis
ibsKing <- snpgdsIBDKING(genofile, sample.id=geno.sample.ids, num.thread=2, type=c("KING-robust"), snp.id=snpset.id0)
ibsKing <- snpgdsIBDSelection(ibsKing)

ggplot(ibsKing[which(ibsKing$IBS0 < 0.1 & ibsKing$kinship > -1),], aes(IBS0, kinship)) + geom_point(alpha=0.7,size=.5) +
  geom_vline(xintercept = 0.001, linetype="dashed") + geom_vline(xintercept = 0.002, linetype="dashed") + geom_vline(xintercept = 0.003, linetype="dashed")
hist(ibsKing$kinship)
# ggplotly()
# //////////////

################
## Further analyses?
# LD decay
# HWE
# Heterozygosity
# //////////////

# BE SURE TO CLOSE THE GDS
snpgdsClose(genofile)
# clean up
rm(genofile, geno.sample.ids, ibsKing, pcaLDP, pc.percentLDP, snpset.id, snpset.id0, snpset)
# //////////////////////////////////////////////////////////



############################################################
## GWASTools analysis

## create phenotype dataframe in special class
scanAnnot = ScanAnnotationDataFrame(data= data.frame(scanID = sampleInfo$ID, 
                                                     AVERAGE_LENGTH = sampleInfo$AVERAGE_LENGTH,
                                                     AVERAGE_LENGTH2 = sampleInfo$AVERAGE_LENGTH2,
                                                     PC1LDP = sampleInfo$PC1LDP,
                                                     PC2LDP = sampleInfo$PC2LDP,
                                                     PC3LDP = sampleInfo$PC3LDP,
                                                     pop = sampleInfo$pop,
                                                     stringsAsFactors = F))

# if you want to make a new GDS file with filtered snps, samples
# snpgdsCreateGenoSet("HDRA-G6-4-RDP1-RDP2-NIAS.gds", "HDRA-G6-4-RDP1-test.gds", snp.id=snpset.id0, sample.id = geno.sample.ids)

# open GDS file
gds <- GdsGenotypeReader("genotype_data/HDRA-G6-4-RDP1-filter2.gds")

# make GenotypeData object, which links phenotype data
genoData <-  GenotypeData(gds, scanAnnot = scanAnnot)
genoData

# get positions to filter on
chrom3start = which(getChromosome(genoData) == __)[1]
chrom6start = which(getChromosome(genoData) == __)[1]

####
## linear regression
res <- assocRegression(genoData,
                       outcome=__,
                       model.type="linear",
                       snpStart = chrom3start, 
                       snpEnd = chrom6start)

# manhattan plot
res$log10pval = -log10(res$Wald.pval)
manhattanPlot(res$Wald.pval, chromosome= res$chr, ylim = c(0, max(res$log10pval) + 1))

# qqplot
qqPlot(pval=res$Wald.pval, truncate=TRUE, main="QQ Plot of Wald Test p-values")
# //

####
## linear regression, pop covariate
res2 <- assocRegression(genoData,
                        outcome="AVERAGE_LENGTH",
                        model.type="linear",
                        covar=__,
                        snpStart = chrom3start, 
                        snpEnd = chrom6start)

# manhattan plot
res2$log10pval = -log10(res2$Wald.pval)
manhattanPlot(res2$Wald.pval, chromosome= res2$chr, ylim = c(0, max(res2$log10pval, na.rm = T) + 1))

# qqplot
qqPlot(pval=res2$Wald.pval, truncate=TRUE, main="QQ Plot of Wald Test p-values")
# //

####
## linear regression, with PCA covariates
res3 <- assocRegression(genoData,
                        outcome="AVERAGE_LENGTH",
                        model.type="linear",
                        block.size = 5000,
                        covar=__,
                        snpStart = chrom3start, 
                        snpEnd = chrom6start)

# manhattan plot
res3$log10pval = -log10(res3$Wald.pval)
manhattanPlot(res3$Wald.pval, chromosome= res3$chr, ylim = c(0, max(res3$log10pval, na.rm = T) + 1))

# qqplot
qqPlot(pval=res3$Wald.pval, truncate=TRUE, main="QQ Plot of Wald Test p-values")
# //

# ////////


##########
## logistic regression
resLogistic <- assocRegression(genoData,
                               outcome=__,
                               model.type="logistic",
                               snpStart = chrom3start, 
                               snpEnd = chrom6start)

# manhattan plot
resLogistic$log10pval = -log10(resLogistic$Wald.pval)
manhattanPlot(resLogistic$Wald.pval, chromosome= resLogistic$chr, ylim = c(0, max(resLogistic$log10pval, na.rm = T) + 1))

# qqplot
qqPlot(pval=resLogistic$Wald.pval, truncate=TRUE, main="QQ Plot of Wald Test p-values")
# //

####
## logistic regression, with pop as covariate
resLogistic2 <- assocRegression(genoData,
                                outcome="AVERAGE_LENGTH2",
                                model.type="logistic",
                                covar=c("pop"),
                                snpStart = chrom3start, 
                                snpEnd = chrom6start)

# manhattan plot
resLogistic2$log10pval = -log10(resLogistic2$Wald.pval)
manhattanPlot(resLogistic2$Wald.pval, chromosome= resLogistic2$chr, ylim = c(0, max(resLogistic2$log10pval, na.rm = T) + 1))

# qqplot
qqPlot(pval=resLogistic2$Wald.pval, truncate=TRUE, main="QQ Plot of Wald Test p-values")
# //

####
## linear regression, with PCA covariates
resLogistic3 <- assocRegression(genoData,
                                outcome="AVERAGE_LENGTH2",
                                model.type="logistic",
                                covar=c("PC1LDP", "PC2LDP", "PC3LDP"),
                                snpStart = chrom3start, 
                                snpEnd = chrom6start)

# manhattan plot
resLogistic3$log10pval = -log10(resLogistic3$Wald.pval)
manhattanPlot(resLogistic3$Wald.pval, chromosome= resLogistic3$chr, ylim = c(0, max(resLogistic3$log10pval, na.rm = T) + 1))

# qqplot
qqPlot(pval=resLogistic3$Wald.pval, truncate=TRUE, main="QQ Plot of Wald Test p-values")
# //

close(gds)
# /////////
save(list=ls(), file="Dataset_out.RData")


# load(file="Dataset_out.RData")
############
## What's next?
# MM testing
# check SNP genotypes in raw data (array, seq) - e.g Do array cluster plots look good?
# search for nearby genes
# //////////