#########################################################################################
## Index Analysis ##
# This code snippet was used to generate intronic analysis using
# the Index package
# https://github.com/charitylaw/Intron-reads/tree/master/Analyses
#########################################################################################

# install.packages("BiocManager")
BiocManager::install("sa-lee/analysis-superintronic")

library(index)
library('analysis-superintronic')

exon <- readRDS(system.file("extdata/exon_dge.Rds", package = "index"))
intron <- readRDS(system.file("extdata/intron_dge.Rds", package = "index"))
group <- readRDS(system.file("extdata/group.Rds", package = "index"))


x <- index_analysis(exon, intron, group)
plot_voom(x)
plot_lcpm_cor(x)
plot_index(x)

#########################################################################################
# Analysis of intron exploration ## https://github.com/charitylaw/Intron-reads/tree/master/Analyses
#########################################################################################

setwd('/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Intron_GTFs')
dir.anno <- '/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Intron_GTFs/annotations'
dir.bam <-'/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Backman19_10.18.2023/totalRNA_New_BAMs'
meta <- read.table(file.path(dir.anno, "metadata.txt"), sep= "\t", header= T)

# 10 Jan 2017 (Last updated 7 Jun 2018)
# Charity Law and Albert Zhang

# Fields to define
exon.anno <- "/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Intron_GTFs/annotations/Exon.txt" # an exon saf file
genebody.anno <- "/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Intron_GTFs/annotations/Genebody.txt" # a genebody saf file
bam <- c('ActD_6hrs_R1_S9_sort.bam','ActD_6hrs_R2_S10_sort.bam','ActD_6hrs_R3_S11_sort.bam','Pol1_8hrs_Aux_R1_S1_sort.bam',
         'Pol1_8hrs_Aux_R2_S2_sort.bam','Pol1_8hrs_Aux_R3_S3_sort.bam', 'Pol1_8hrs_Aux_R4_S4_sort.bam','Pol2_8hrs_Aux_R1_S5_sort.bam','Pol2_8hrs_Aux_R2_S6_sort.bam','Pol2_8hrs_Aux_R3_S7_sort.bam',
         'Pol2_8hrs_Aux_R4_S8_sort.bam','WT_HCT116_CTRL_R1_S12_sort.bam','WT_HCT116_CTRL_R2_S13_sort.bam','WT_HCT116_CTRL_R3_S14_sort.bam') # bam file names (excluding 1 pol1 and 1 pol2 rep)
is.paired.end <- FALSE # seq protocol for read summarisation
sample.names <- meta$SampleNames # sample names (unique)
groups <- meta$Group # groups for samples (non-unique)
threads = 12

# Load libraries
library(limma)
library(edgeR)
library(Rsubread)
library(stringr)

# Create folders
dir.create("Counts", showWarnings=FALSE)
dir.create("Results", showWarnings=FALSE)
dir.create("Plots", showWarnings=FALSE)


### SUMMARISE READS TO COUNTS

# Exon
anno <- read.delim(exon.anno)
counts <- featureCounts(annot.ext=anno, files=paste0(dir.bam, '/',bam), isPairedEnd=is.paired.end, useMetaFeatures=TRUE, allowMultiOverlap=FALSE, strandSpecific=0, nthreads=threads)
colnames(counts$counts) <- colnames(counts$stat)[-1] <- sample.names
save(counts, file="Counts/Exon.RData")  

# Genebody
anno <- read.delim(genebody.anno)
counts <- featureCounts(annot.ext=anno, files=paste0(dir.bam, '/',bam), isPairedEnd=is.paired.end, useMetaFeatures=TRUE, allowMultiOverlap=FALSE, strandSpecific=0, nthreads=threads)
colnames(counts$counts) <- colnames(counts$stat)[-1] <- sample.names
save(counts, file="Counts/Genebody.RData")    

# Intron  
load("Counts/Exon.RData")
exon <- counts
load("Counts/Genebody.RData")
genebody <- counts
counts$counts <- pmax(genebody$counts-exon$counts,0)
save(counts, file="Counts/Intron.RData")


### CALC READ PERCENTAGES

# Get library size and number of exon counts
load("Counts/Exon.RData")
ntotal <- colSums(counts$stat[,-1])
nexon <- colSums(counts$counts)

# Get total intron reads
load("Counts/Intron.RData")
nintron <- colSums(counts$counts) 

# Calculate read percentages
read.percentages <- rbind(Exon=nexon, Intron=nintron)/rep(ntotal, each=2)
colnames(read.percentages) <- groups

# Save
write.table(read.percentages, file="Results/read_percentage.txt", sep="\t")

# Plot
group.f <- as.factor(groups)
pdf("Plots/read_percentage.pdf")
par(mar=c(10,5,1,1))
plot(as.numeric(group.f)+rnorm(length(group.f), sd=0.05), read.percentages["Exon",],
     ylim=c(0,1), col="dodgerblue", ylab="Proportion of reads", xlab="", xaxt="n")
points(as.numeric(group.f)+rnorm(length(group.f), sd=0.05), read.percentages["Intron",], col="magenta")
axis(side=1, at=1:nlevels(group.f), labels=levels(group.f), las=2)
legend("topleft", fill=c("dodgerblue", "magenta"), legend=rownames(read.percentages), bty="n")
dev.off()


# MDS PLOTS

# Create and save MDS plots
pdf("Plots/sample_clusters.pdf", width=8, height=4)
par(mar=c(5,5,2,1))
par(mfrow=c(1,2))
# Make exon MDS plot
load("Counts/Exon.RData")
lcpm <- cpm(counts$counts, log=TRUE)
mds.exon <- plotMDS(lcpm, labels=groups, col=as.numeric(as.factor(groups)), main="Exon")

# Make intron MDS plot
load("Counts/Intron.RData")
lcpm <- cpm(counts$counts, log=TRUE)
mds.intron <- plotMDS(lcpm, labels=groups, col=as.numeric(as.factor(groups)), main="Intron")
dev.off()
# Save 
save(mds.exon, mds.intron, file="Results/sample_clusters.RData")


# GENE SIGNAL

# Exon counts
load("Counts/Exon.RData")
exon <- DGEList(counts$counts, genes=counts$annotation)
# Intron counts
load("Counts/Intron.RData")
intron <- DGEList(counts$counts, genes=counts$annotation)
m <- match(exon$genes$GeneID, intron$genes$GeneID)
intron <- intron[m,] 

# Genes with signal in exon and intron regions
exon.signal <- exon$counts>=3
intron.signal <- intron$counts>=3
signal <- matrix(NA, nrow=nrow(exon), ncol=ncol(exon))
signal[exon.signal & intron.signal] <- "Exon and intron signal"
signal[!exon.signal & !intron.signal] <- "No signal"
signal[exon.signal & !intron.signal] <- "Exon signal"
signal[!exon.signal & intron.signal] <- "Intron signal"

# Mark genes with a single exon
nexons <- str_count(exon$genes$Chr, ";")+1
single.exon <- rep("", length(nexons))
single.exon[nexons==1] <- " (no introns)"
signal <- matrix(paste0(signal, single.exon), nrow=length(single.exon))
colnames(signal) <- sample.names
rownames(signal) <- exon$genes$GeneID

# Summarise 
signal.summarised <- apply(signal, 2, function(x) table(x))
signal.summarised <- signal.summarised/nrow(signal)

# Save
write.table(signal, file="Results/gene_signal.txt", sep="\t")
write.table(signal.summarised, file="Results/gene_signal_summarised.txt", sep="\t")

# Plot
dir.create("Plots", showWarnings=FALSE)
group <- as.factor(groups)
pdf("Plots/gene_signal.pdf")
par(mar=c(10,5,1,1))
MAX <- max(signal.summarised)+0.01
plot(as.numeric(group)+rnorm(length(group), sd=0.05), signal.summarised["Exon and intron signal",],
     ylim=c(0,MAX), col="grey", ylab="Proportion of genes", xlab="", xaxt="n")
points(as.numeric(group)+rnorm(length(group), sd=0.05), signal.summarised["Exon signal",], col="dodgerblue")
points(as.numeric(group)+rnorm(length(group), sd=0.05), signal.summarised["Intron signal",], col="magenta")
points(as.numeric(group)+rnorm(length(group), sd=0.05), signal.summarised["Exon signal (no introns)",], col="limegreen")
axis(side=1, at=1:nlevels(group), labels=levels(group), las=2)
legend("topleft", fill=c("grey", "dodgerblue", "limegreen", "magenta"), legend=c("Exon and intron signal", "Exon signal", "Exon signal (no introns)", "Intron signal"), bty="n")
dev.off()


# EXON VS INTRON LOG-COUNTS

load("Counts/Exon.RData")
exon <- DGEList(counts$counts, genes=counts$annotation)
# Intron counts
load("Counts/Intron.RData")
intron <- DGEList(counts$counts, genes=counts$annotation)
m <- match(exon$genes$GeneID, intron$genes$GeneID)
intron <- intron[m,] 
signal <- read.delim("Results/gene_signal.txt", header=TRUE)
logcounts.signal <- rep(NA, ncol(exon))
names(logcounts.signal) <- groups
for (i in 1:length(sample.names)){
  exon.i <- log2(exon$counts[,i]+1)
  intron.i <- log2(intron$counts[,i]+1)
  signal.i <- signal[,i]=="Exon and intron signal"
  logcounts.signal[i] <- cor(exon.i[signal.i], intron.i[signal.i])
}
pdf("Plots/gene_signal_correlation.pdf", height=4, width=4)
par(mar=c(10,4,1,1))
boxplot(logcounts.signal~groups, las=2)
dev.off()
# Save
write.table(logcounts.signal, file="Results/gene_signal_correlation.txt", row.names=FALSE)

#########################################################################################
# Index analysis 
#########################################################################################

## directories
setwd('/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Intron_GTFs')
dir.anno <- '/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Intron_GTFs/annotations'
dir.bam <-'/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Backman19_10.18.2023/totalRNA_New_BAMs'

library(index)
library(magrittr)
library(stringr)

## get meta data 
meta <- read.table(file.path(dir.anno, "metadata.txt"), sep= "\t", header= T)
sample.names <- meta$SampleNames # sample names (unique)
group <- meta$Group # groups for samples (non-unique)

##-------------------------------------------------------## DE analysis
## Helper function to get contrasts for each comparitor
get_pairwise_contrasts <- function(design, cols = ncol(design)) {
  group_names <- colnames(design)
  
  contrasts <- list()
  contrast_names <- character()
  for (i in 1:(cols - 1)) {
    for (j in (i + 1):cols) {
      contr <- numeric(length = cols)
      contr[i] <- 1
      contr[j] <- -1
      contrast_name <- glue::glue("{group_names[i]} vs {group_names[j]}")
      contrasts <- append(contrasts, list(contr))
      contrast_names <- append(contrast_names, contrast_name)
    }
  }
  setNames(contrasts, contrast_names)
}

## Function to get results using INDEX for all comparisons
getDGE <- function(vec, group){
  
  ## this function only works for single comparisons
  ## vec is the order of sample columns in the counts data. must match group list
  
  ## Load data
  exon <- load("Counts/Exon.RData")
  exon <- DGEList(counts$counts[,vec], genes=counts$annotation)
  intron <- load("Counts/Intron.RData")
  intron <- DGEList(counts$counts[,vec], genes=counts$annotation)
  genebody <- load("Counts/Genebody.RData")
  genebody <- DGEList(counts$count[,vec], genes=counts$annotation)
  group <- group
  
  ## get intron lengths
  intron$genes$Length <- genebody$genes$Length - exon$genes$Length + 1
  
  ## get experiment design
  design <- model.matrix(~ 0 + group) %>%
    set_colnames(colnames(.) %>% str_remove("group"))
  print(paste0("design: "))
  print(design)
  
  ## get pairwise contrast for design
  pairwise_contrasts <- get_pairwise_contrasts(design)
  print(paste0("pairwise contrasts: "))
  print(pairwise_contrasts)
  
  ## use INDEX to get intron diff expression
  res <-index_analysis(exon, intron, group = group, design = design, contrast = pairwise_contrasts[[1]], p.value = 0.05)
  return(res)
  
}

actd.res <- getDGE(vec = c(12,13,14,1,2,3), group = c("CTRL_0hrs", "CTRL_0hrs", "CTRL_0hrs","ActD_6hrs", "ActD_6hrs", "ActD_6hrs"))
pol1.res <- getDGE(vec = c(12,13,14,4,5,6,7), group = c("CTRL_0hrs", "CTRL_0hrs", "CTRL_0hrs","Pol1_6hrs", "Pol1_6hrs", "Pol1_6hrs", "Pol1_6hrs"))
pol2.res <- getDGE(vec = c(12,13,14,8,9,10,11), group = c("CTRL_0hrs", "CTRL_0hrs", "CTRL_0hrs","Pol2_6hrs", "Pol2_6hrs", "Pol2_6hrs", "Pol2_6hrs"))

# Save
write.table(pol1.res$tops$intron, file="Results/pol1.intron.dge.txt", row.names=FALSE)
write.table(pol2.res$tops$intron, file="Results/pol2.intron.dge.txt", row.names=FALSE)
write.table(actd.res$tops$intron, file="Results/actd.intron.dge.txt", row.names=FALSE)
