
# upregulated color: "#D46780"
# downregulated color: "#798234"

rm(list = ls())
gc()

## Directories
## Set your root directory here
root<- "/Volumes/external hd/IBiS/Backman_Lab/Transcription Publication/repositories/dryad/"
dir.rsem.txi <- "/RNAseq/counts/RSEM/isoforms"
dir.rsem <- "/RNAseq/counts/RSEM"
dir.htseq <- "/RNAseq/counts/HTseq"

# import packages
require("DESeq2")
require("ggplot2")
require("dplyr")
require("magrittr")
require("tximport")
require("forcats")
require("stringr")
require("matrixStats")
require("ggrepel")
require("pheatmap")
require("viridis")
require("dendsort")
require("tidyr")
require("ggpubr")
require("ggside")

#########################################################################################
# Functions
#########################################################################################

## A function to filter results for top DE genes
Filter_Results <- function(results, lfc.up=0.58, lfc.dwn=-0.58,pvadj=0.05) {
  
  resultsDN <- as.data.frame(results[which(results$log2FoldChange < lfc.dwn),])
  resultsDN <- as.data.frame(resultsDN[which(resultsDN$padj <= pvadj),])
  length(rownames(resultsDN)) ## Gives count of down regulated
  
  resultsUP <- as.data.frame(results[which(results$log2FoldChange  > lfc.up),])
  resultsUP <- as.data.frame(resultsUP[which(resultsUP$padj <= pvadj),])
  length(rownames(resultsUP)) ## Gives count of up regulated  
  
  ## add a column of NAs
  results$diffexpressed <- "NO"
  ## If log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  results$diffexpressed[results$log2FoldChange > lfc.up & results$padj <= pvadj] <- "Upregulated"
  ## If log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  results$diffexpressed[results$log2FoldChange  < lfc.dwn & results$padj <= pvadj] <- "Downregulated"
  
  ## Make counts vectors
  results$up <- ifelse(results$diffexpressed == "Upregulated", 1, 0)
  results$down <- ifelse(results$diffexpressed == "Downregulated", 1, 0)
  
  ## Make lfc vectors
  results$up.lfc <- ifelse(results$diffexpressed == "Upregulated", results$log2FoldChange, 0)
  results$down.lfc <- ifelse(results$diffexpressed == "Downregulated", results$log2FoldChange, 0)
  
  ## Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  results$delabel <- NA
  results$delabel[results$diffexpressed != "NO"] <- results$symbol[results$diffexpressed != "NO"]
  
  return(results)
  
}


## A function to annotate DGE results. If using transcript level quantification instead of gene level quant, set transcript=TRUE
Annotate_results <- function(results, filename, transcript,
                             dir.res = "RNAseq/results/DGE_counts"
){
  
  require("AnnotationHub")
  require("AnnotationDbi")
  
  ah <- AnnotationHub() ## Connect to Annotation Hub
  
  db <- query(ah, c("EnsDb","Homo sapiens","111"))[[1]] ## Make a query to your organism/database
  results <- data.frame(results)
  
  if (transcript == F) {
    
    annotations <- genes(db, return.type = "data.frame")
    annot <- annotations %>%
      dplyr::select(gene_id, canonical_transcript, symbol,gene_biotype, seq_name, gene_seq_start, gene_seq_end) %>%
      dplyr::filter(gene_id %in% rownames(results)) ## genes.results
    
    ## Order results to match them
    order <- match(rownames(results), annot$gene_id) ## genes.results
    annot <- annot[order,]
    
    ## Annotations
    results$gene_id <- annot$gene_id
    results$transcript <- annot$canonical_transcript
    results$symbol <- annot$symbol
    results$gene_biotype <- annot$gene_biotype
    results$seq_name <- annot$seq_name
    results$gene_seq_start <- annot$gene_seq_start
    results$gene_seq_end <- annot$gene_seq_end
    
  } else {
    
    keys <- keys(db, c("TXID"))
    annotations <- AnnotationDbi::select(db, keys, c("GENEID", "TXNAME", "SYMBOL", "GENENAME","GENEBIOTYPE","TXBIOTYPE", "SEQNAME", "TXSEQSTART","TXSEQEND","EXONID","EXONSEQSTART", "EXONSEQEND"), "TXID")
    
    annot <- annotations %>%
      dplyr::select("GENEID", "TXID", "TXNAME", "SYMBOL", "GENENAME","GENEBIOTYPE","TXBIOTYPE", "SEQNAME", "TXSEQSTART","TXSEQEND","EXONID","EXONSEQSTART", "EXONSEQEND") %>%
      dplyr::filter(TXID %in% rownames(results)) ## genes.results
    
    ## Order results to match them
    order <- match(rownames(results), annot$TXID) ## genes.results
    annot <- annot[order,]
    
    ## Annotations
    results$gene_id <- annot$GENEID
    results$transcript <- annot$TXID
    results$symbol <- annot$GENENAME
    results$gene_biotype <- annot$GENEBIOTYPE
    results$tx_biotype <- annot$TXBIOTYPE
    results$seq_name <- annot$SEQNAME
    results$tx_seq_start <- annot$TXSEQSTART
    results$tx_seq_end <- annot$TXSEQEND
    results$exon_id <- annot$EXONID
    results$exon_seq_start <- annot$EXONSEQSTART
    results$exon_seq_end <- annot$EXONSEQEND
    
  }
  
  results$padj <- ifelse(is.na(results$padj), 1, results$padj) ## NA filtering, padj
  results$pvalue <- ifelse(is.na(results$pvalue), 1, results$pvalue) ## NA filtering, pv
  results <- na.omit(results)
  
  ## Order results and save Files as CSV function
  results <- results[order(results$padj),]
  results$Count <-rep(1, nrow(results))
  
  write.csv(as.data.frame(results), file=file.path(root,dir.res,paste0(filename, ".csv")))
  
  return(results)
  
}

## This function is for annotating by genomic features like promoters, exons, etc
## mapping geneIDs blog: https://medium.com/computational-biology/gene-id-mapping-using-r-14ff50eec9ba

## A function for annotating genes with additional information
annotate.peak <- function(peaks.file){
  
  require(ChIPseeker)
  require(TxDb.Hsapiens.UCSC.hg38.refGene)
  require(plyranges)
  
  res <- as.data.frame(peaks.file)
  res.gr<- res %>% 
    #transform( seqnames= paste0("chr",res$seq_name), start = res$gene_seq_start, end = res$gene_seq_end)  %>% 
    transform( seqnames= paste0("chr",res$seq_name), start = res$tx_seq_start, end = res$tx_seq_end)  %>% 
    as_granges()
  
  
  ## Prepare annotation objects for mapping
  txdb <- TxDb.Hsapiens.UCSC.hg38.refGene ## Refgene and knowngene dbs have signficant differences
  ## see: https://vatlab.github.io/vat-docs/applications/annotation/genes/refgene/
  
  ## Call annotation
  anno <- annotatePeak( res.gr, TxDb = txdb , tssRegion=c(-3000, 3000))
  plotAnnoPie(anno)
  anno <- data.frame(anno@anno)
  
  return(anno)
  
}

## A function that calls on a function from ggpubr to generate an MA plot from DESeq results object
plot_MA <- function(coef, title, results, dncol="#D46780", upcol="#798234"){
  
  require("apeglm")
  require("ggpubr")
  
  # shrink noisy log2 fold change estimates for coefficients: DESeq2 function
  res <- lfcShrink(dds, coef=coef, type="apeglm")
  res <- res[rownames(res)%in% rownames(results),]
  
  g1<-ggmaplot(res, 
               main = title,
               fdr = 0.10, fc = 1, size = 2,
               palette = c(upcol, dncol, "darkgray"),
               legend = "top", top = 0,
               alpha=0.8,
               xlab = "Log2 mean expression",
               ylab = "Log2 fold change",
               font.legend = c(16, "black", "Arial"),
               font.main = c(16,"black", "Arial"),
               ggtheme = theme_classic()
  ) + theme(text = element_text(size = 16, color= "black", family = "Arial"))
  
  return(g1)
  
}

## Function to generate volcano plot
DE_Vol_Plot <- function(results, plt.title, lfc.up=0.58, lfc.dwn=-0.58,pvadj=0.05,
                        y.limits, x.limits=c(-7,7),dncol="#D46780", upcol="#798234") {
  require(ggrepel)
  require(ggplot2)
  
  ## Generate data frame
  genes.rna <- data.frame(symbol = results$symbol,lfc = results$log2FoldChange, padj=results$padj, base=results$baseMean)
  
  resultsDN <- as.data.frame(genes.rna[which(genes.rna$lfc < lfc.dwn),])
  resultsDN <- as.data.frame(resultsDN[which(resultsDN$padj < pvadj),])
  length(rownames(resultsDN)) ## Gives count of down regulated
  
  resultsUP <- as.data.frame(genes.rna[which(genes.rna$lfc > lfc.up),])
  resultsUP <- as.data.frame(resultsUP[which(resultsUP$padj < pvadj),])
  length(rownames(resultsUP)) ## Gives count of up regulated  
  
  ## add a column of NAs
  genes.rna$diffexpressed <- "NO"
  ## If log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  genes.rna$diffexpressed[genes.rna$lfc > lfc.up & genes.rna$padj < pvadj] <- "Upregulated"
  ## If log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  genes.rna$diffexpressed[genes.rna$lfc < lfc.dwn & genes.rna$padj < pvadj] <- "Downregulated"
  
  ## Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  genes.rna$delabel <- NA
  genes.rna$delabel[genes.rna$diffexpressed != "NO"] <- genes.rna$symbol[genes.rna$diffexpressed != "NO"]
  
  ## Jitter labeling
  pos <- position_jitter(width = 0.2, seed = 1)
  
  vplot <- genes.rna %>% mutate(mean.intensity = base/max(base)) %>% 
    ggplot(aes(x=lfc, y=-log10(padj), col=diffexpressed, size = mean.intensity #, label=delabel
    )) +
    geom_point(position = pos, alpha=0.7)+ theme_classic()+
    scale_color_manual(name="Upregulated", values=c(dncol, "lightgrey", upcol))+
    geom_vline(xintercept=c(lfc.dwn, lfc.up), col="red", lty = "dashed") +
    geom_hline(yintercept=-log10(pvadj), col="red", lty = "dashed")+ 
    #geom_text_repel(segment.color = 'transparent', position = pos, size = 3, xlim  = c(-6,6)) +
    labs(title=plt.title, subtitle=paste(paste("Downregulated:",length(rownames(resultsDN))), paste("| Upregulated:",length(rownames(resultsUP))))) +
    scale_y_continuous(name="-Log 10 Padj", limits=c(0,y.limits), breaks=seq(0,y.limits, by =50))+
    scale_x_continuous(name="Log 2 Fold Change", limits=x.limits, breaks=c(-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7))+
    theme( legend.position="none",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) 
  
  return(vplot)
  
}
#########################################################################################
# Set Directory, pass in  counts, filter counts for DGE analysis (RSEM counts)
#########################################################################################

## Get sample metadata
sample <- read.table(file.path(root,dir.rsem, "metadata.csv"), sep= ",", header= TRUE)
newcolnames <-sample$Alias
sample

## Import files
files <- file.path(root,dir.rsem, paste0(sample$ID))
files

## Imports RSEM TPM files into a list genes.results
cts <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

##-------------------------------------------------------## Get transcript level isoforms.results

files <- file.path(root,dir.rsem.txi, paste0(str_replace_all(sample$ID, "genes", "isoforms")))
files

cts<- tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE)
head(cts)

##-------------------------------------------------------## write out counts table for GEO submission

colnames(cts$counts) <- c("ActD_Rep1", "ActD_Rep2", "ActD_Rep3", "Pol1_Rep1", "Pol1_Rep2", "Pol1_Rep3", 
  "Pol2_Rep1", "Pol2_Rep2", "Pol2_Rep3", "ctrl_Rep1", "ctrl_Rep2", "ctrl_Rep3")

dir.res = "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/results/DGE_counts"
write.table(as.data.frame(cts$counts), file=paste0(root,dir.res,"/", "rawcounts.tsv"), sep = "\t")

##-------------------------------------------------------## 

## Move col names from sample info to counts matrix
colnames(cts$counts)<-newcolnames

#Replace any transcripts with length 0 with a 1
cts$length[cts$length == 0] <- 1
head(cts$counts)

# logical test to see if ^ is true now
identical(sample$Alias, colnames(cts$counts))

# Store the file names in the sample info matrix as row names
row.names(sample) <- sample$Alias
head(sample, n=10)

#factor out dependencies
sample$Conditions <- factor(sample$Conditions, levels = c("CTRL_0hrs","ActD_6hrs","Pol1_6hrs","Pol2_6hrs"))
sample<- mutate(sample, Conditions = fct_relevel(Conditions, c( "CTRL_0hrs","ActD_6hrs","Pol1_6hrs","Pol2_6hrs")))
sample$Label <- factor(sample$Label, levels = dput(as.character(unique(sample$Label))))

#########################################################################################
# Set Directory, pass in  counts, filter counts for DGE analysis (HTSEQ counts)
#########################################################################################

list.files(file.path(root,dir.htseq))

## Get counts
cts <- read.table(file.path(root,dir.htseq,"htseq_all_sample_count.tsv"),sep="\t",header= TRUE, stringsAsFactors = FALSE)
head(cts)

## remove tallies and rownames column
cts <- cts[-c(1:5),]
rownames(cts) <- cts$X
cts <- cts[,-c(1)]

## get sample metadata
sample <- read.table(file.path(root, dir.rsem, "metadata.csv"), sep= ",", header= TRUE)
names <- colnames(cts)
sample$ID2 <- names

newcolnames <-sample$Alias

# move col names from sample info to counts matrix
colnames(cts)<-newcolnames
colnames(cts)

# match the sample information and count matrix together
cts <- cts[,match(sample$Alias, colnames(cts))]
# logical test to see if ^ is true now
identical(sample$Alias, colnames(cts))

cts<- data.frame(as.matrix(cts))

# Store the file names in the sample info matrix as row names
row.names(sample) <- sample$Alias
head(sample, n=10)

#factor out dependencies
sample$Conditions <- factor(sample$Conditions, levels = c( "CTRL_0hrs","ActD_6hrs","Pol1_6hrs","Pol2_6hrs"))
sample$Label <- factor(sample$Label, levels = dput(as.character(unique(sample$Label))))


#########################################################################################
# Create DESeq data object and pply transformations to data for pre-DGE QC
#########################################################################################

## creat DDS object using HTseq counts
dds <- DESeqDataSetFromMatrix(countData = cts, colData = sample, design = ~ Conditions)

## Create DDS object using RSEM counts
dds <- DESeqDataSetFromTximport(cts, colData = sample, design = ~ Conditions)

## check dds object
summary(dds)

## Pre-filtering the dataset
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 0, ]
nrow(dds)

## Variance stabilizing transformation (VST) | regularized log transformation | and normalized counts transformation to make data more homoskedastic
vst <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
ntd <-normTransform(dds, f=log2, pc=1)

## Add another column which defines the group that each sample is in
dds$group <- factor(paste0(dds$Conditions))

## Find size factors to normalize differences in sequencing depth across samples
dds <- estimateSizeFactors(dds)

#########################################################################################
# Base R PCA Plots
#########################################################################################

## Count table for PCA
cts.norm <- as.matrix(round(counts(dds)))
head(cts.norm)

scaled = DESeq2::rlog(cts.norm, blind = TRUE)

## Next we calculate the per-row variance in the transformed table.
vars = rowVars(scaled)

## Then we re-order the rows so that the most variable rows are on top.
scaled = scaled[order(vars, decreasing = TRUE),]

## ## Here we perform PCA using `prcomp`. We transpose (`t`) the matrix and only operate on the entire data set
pca = prcomp(t(scaled))

## We calculate the percent variance explained by each principal component by extracting these values from the summary.
percentVar = round(100 * summary(pca)$importance[2,])
cumuVar = round(100 * summary(pca)$importance[3,])

## Look at top variance
head(pca$x, n=12)

prTab = data.frame(pca$x, colData(dds)) ## All conditions

prTab %>% #PCA1_2_AllSamples_withlabels ## pca.pol1only ## pca.pol1only.nolabels ## pca.1.2.nopol2 ## pca.1.2.nopol2.nolabels
  ggplot(aes(PC1, PC2, color = Conditions)) +
  geom_point(size = 8, alpha = 1) + theme_classic()+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  scale_color_manual(values = c("black", "darkblue","slateblue", "steelblue3")) +
  coord_fixed() + geom_text_repel(segment.color = 'transparent', aes(label = NA), size = 3, family = "Arial Narrow", hjust = 0)+
  theme(legend.position="bottom", plot.title = element_text(size=20), text = element_text(size=20, family="Arial"))

## Plot top PCAs for PCA2 and 3 
prTab %>% #PCA2_3_AllSamples_NOlabels ## pca.2.3.pol1only ## pca.2.3.pol1only.nolabels ## pca.2.3.nopol2 ## pca.2.3.nopol2.nolabels
  ggplot(aes(PC2, PC3, color = Conditions)) +
  geom_point(size = 5, alpha = 1) + theme_classic()+
  xlab(paste0("PC2: ", percentVar[2], "% variance")) + ylab(paste0("PC3: ", percentVar[3], "% variance")) + 
  scale_color_manual(values = c("black", "darkblue","slateblue", "steelblue3")) +
  coord_fixed() + geom_text_repel(aes(label = NA), size = 3, family = "Arial Narrow", hjust = 0)+
  theme(legend.position="bottom", plot.title = element_text(size=20), text = element_text(size=20, family="Arial"))
  
## Plot cumulative variance
b = barplot(percentVar, horiz = FALSE, las = 1, ylab = 'Variance (%)', ylim = c(0,105))
lines(b, cumuVar, type = 'b', col = 'red')

## Plot Eigens
eigenvalues = pca$sdev^2
plot(c(1:12), eigenvalues, xlab = 'Principal Components', ylab = 'eigenvalue', las = 1)

## determing which PCAs to keep based on: 
which(eigenvalues > mean(eigenvalues))

#########################################################################################
# Figure 3H: Generate Heatmap for QC of count data
#########################################################################################

## Hierachal clustering map plot for sample QC
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

## Extract the rlog matrix from the object
rld_mat <- assay(rld) ## All samples 

## Compute pairwise correlation values
rld_cor <- cor(rld_mat) 

colnames(rld_cor) <- c("ActD R1", "ActD R2", "ActD R3", "Pol1 R1", "Pol1 R2", "Pol1 R3", "Pol2 R1", "Pol2 R2", "Pol2 R3", "ctrl R1", "ctrl R2", "ctrl R3")
rownames(rld_cor) <- c("ActD R1", "ActD R2", "ActD R3", "Pol1 R1", "Pol1 R2", "Pol1 R3", "Pol2 R1", "Pol2 R2", "Pol2 R3", "ctrl R1", "ctrl R2", "ctrl R3")

## organize hierachal clusters
mat <- rld_cor
mat_cluster_cols <- hclust(dist(t(mat)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(mat)))
length(mat)

## Color scales
colors <- colorRampPalette(c("#798234","white","#D46780"))(length(mat))

# plot heat map ## sample.corr.plot.nopol2 ## sample.corr.plot.pol1only
p<-pheatmap(mat, 
            color = colors,
            fontsize = 16,
            fontsize_row = 16, 
            fontsize_col = 16,
            cellwidth = 20,
            cellheight = 20,
            cluster_cols= mat_cluster_cols,
            cluster_rows= mat_cluster_rows,
            main = "",
            annotation_legend = TRUE,
            scale="row")
p

#########################################################################################
# Differential expression analysis and 'results' matrices
#########################################################################################

## Generate DGE DESeqDataSet data frame
dds <- DESeq(dds)

## Which contrasts to use?
resultsNames(dds)

## Generate a results matrix for Treatment vs. Control with an alpha of 0.05
res.actd <- results(dds, contrast = c("Conditions","ActD_6hrs", "CTRL_0hrs"), alpha = 0.05)
res.pol1 <- results(dds, contrast = c("Conditions", "Pol1_6hrs", "CTRL_0hrs"), alpha = 0.05)
res.pol2 <- results(dds, contrast = c("Conditions", "Pol2_6hrs", "CTRL_0hrs"), alpha = 0.05)

#########################################################################################
# Annotate the results using Annotation Hub and read results out
#########################################################################################

## Annotate genes ## When using gene RSEM counts
res.actd <- Annotate_results(res.actd, "Conditions_ActD_6hrs_vs_CTRL_0hrs.genes", transcript=F)
res.pol1 <- Annotate_results(res.pol1, "Conditions_Pol1_6hrs_vs_CTRL_0hrs.genes", transcript=F)
res.pol2 <- Annotate_results(res.pol2, "Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes", transcript=F)

## Annotate transcripts ## When using transcript RSEM counts
res.actd <- Annotate_results(res.actd, "Conditions_ActD_6hrs_vs_CTRL_0hrs.tx", transcript=T)
res.pol1 <- Annotate_results(res.pol1, "Conditions_Pol1_6hrs_vs_CTRL_0hrs.tx", transcript=T)
res.pol2 <- Annotate_results(res.pol2, "Conditions_Pol2_6hrs_vs_CTRL_0hrs.tx", transcript=T)

#####################################################################################
# Figure 3D-F: Get counts of transcript biotype and cum lfc                            
#####################################################################################

## generate table of up/down values
res.pol1.filtered <- Filter_Results(res.pol1)
res.pol2.filtered <- Filter_Results(res.pol2)
res.actd.filtered <- Filter_Results(res.actd)

##-------------------------------------------------------## transcript Annotation

pol1.anno <- annotate.peak(res.pol1.filtered)
pol2.anno <- annotate.peak(res.pol1.filtered)
actd.anno <- annotate.peak(res.actd.filtered)

##-------------------------------------------------------## Plot Figures 3D-F

## get table of number of DEGs up and down
feat_tbl.1 <-res.pol1.filtered  %>% group_by(tx_biotype) %>% 
  reframe(down = as.numeric(-sum(down)), up = as.numeric(sum(up))) 
  
feat_tbl.2 <- res.actd.filtered  %>% group_by(tx_biotype) %>% 
  reframe(down = as.numeric(-sum(down)), up = as.numeric(sum(up))) 

feat_tbl.3 <- res.pol2.filtered  %>% group_by(tx_biotype) %>% 
  reframe(down = as.numeric(-sum(down)), up = as.numeric(sum(up))) 

## clean up data for plotting
feat_tbl.1$name <- rep("Pol 1", nrow(feat_tbl.1))
feat_tbl.2$name <- rep("ActD", nrow(feat_tbl.2))
feat_tbl.3$name <- rep("Pol 2", nrow(feat_tbl.3))

## Gather Pol 1 data for plotting
tbl <- gather(feat_tbl.1, diff, sum, 2:3)
tbl <- tbl[tbl$sum != 0,]
tbl <- tbl[abs(tbl$sum) > 25,]

## # deg.summary.pol1 ## deg.tx.summary.pol1
tbl %>% 
  ggplot(aes(x=reorder(tx_biotype, -abs(sum)), y=as.numeric(sum) , fill=diff, color = "black"), group_by(tx_biotype)) +
  scale_y_continuous(breaks = c(seq(min(as.numeric(tbl$sum)), 0, by =1000), seq(0, max(as.numeric(tbl$sum)),by =1000))  ,limits=c(min(as.numeric(tbl$sum)), max(as.numeric(tbl$sum)))) + 
  geom_bar(stat="identity", width=0.75, position = "stack", alpha = 1) +
  scale_color_manual(values=c("black","black"))+
  scale_fill_manual(values=c("#D46780","#798234")) + theme_classic()+
  scale_x_discrete(expand = c(0.09, 0.09))+
  theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position="", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"))+
  guides(fill = guide_legend(title = "Transcript Type"))+ ylab("Number of DEGs") + xlab("")+ ggtitle(paste0("Pol 1"))

## Gather Pol 2 data for plotting
tbl <- gather(feat_tbl.2, diff, sum, 2:3)
tbl <- tbl[tbl$sum != 0,]
tbl <- tbl[abs(tbl$sum) > 25,]

## # deg.summary.pol2 ## deg.tx.summary.pol2
tbl %>% 
  ggplot(aes(x=reorder(tx_biotype, -abs(sum)), y=as.numeric(sum) , fill=diff, color = "black"), group_by(tx_biotype)) +
  scale_y_continuous(breaks = c(seq(min(as.numeric(tbl$sum)), 0, by =1000), seq(0, max(as.numeric(tbl$sum)),by =1000))  ,limits=c(min(as.numeric(tbl$sum)), max(as.numeric(tbl$sum)))) + 
  geom_bar(stat="identity", width=0.75, position = "stack", alpha = 1) +
  scale_color_manual(values=c("black","black"))+
  scale_fill_manual(values=c("#D46780","#798234")) + theme_classic()+
  scale_x_discrete(expand = c(0.09, 0.09))+
  theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position="", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"))+
  guides(fill = guide_legend(title = "Transcript Type"))+ ylab("Number of DEGs") + xlab("")+ ggtitle(paste0("Pol 2"))

## Gather ActD data for plotting
tbl <- gather(feat_tbl.3, diff, sum, 2:3)
tbl <- tbl[tbl$sum != 0,]
tbl <- tbl[abs(tbl$sum) > 25,]

## # deg.summary.actd ## deg.tx.summary.actd
tbl %>% 
  ggplot(aes(x=reorder(tx_biotype, -abs(sum)), y=as.numeric(sum) , fill=diff, color = "black"), group_by(tx_biotype)) +
  scale_y_continuous(breaks = c(seq(min(as.numeric(tbl$sum)), 0, by =1000), seq(0, max(as.numeric(tbl$sum)),by =1000))  ,limits=c(min(as.numeric(tbl$sum)), max(as.numeric(tbl$sum)))) + 
  geom_bar(stat="identity", width=0.75, position = "stack", alpha = 1) +
  scale_color_manual(values=c("black","black"))+
  scale_fill_manual(values=c("#D46780","#798234")) + theme_classic()+
  scale_x_discrete(expand = c(0.09, 0.09))+
  theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position="", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"))+
  guides(fill = guide_legend(title = "Transcript Type"))+ ylab("Number of DEGs") + xlab("")+ ggtitle(paste0("ActD"))

#####################################################################################
# Figure 3G: Build scatter plots of log fold change and counts                                     
#####################################################################################

##-------------------------------------------------------## For transcripts
res.pol1.filtered$tx_id <- rownames(res.pol1.filtered)
res.pol2.filtered$tx_id <- rownames(res.pol2.filtered)
res.actd.filtered$tx_id <- rownames(res.actd.filtered)

pol1.df <- dplyr::select(res.pol1.filtered, tx_id, log2FoldChange, padj, gene_biotype, diffexpressed)
pol2.df <- dplyr::select(res.pol2.filtered, tx_id, log2FoldChange, padj, gene_biotype, diffexpressed)
actd.df <- dplyr::select(res.actd.filtered, tx_id, log2FoldChange, padj, gene_biotype, diffexpressed)

plot.1 <- merge(pol1.df, pol2.df, by = "tx_id")
plot.2 <- merge(pol1.df, actd.df, by = "tx_id")
plot.3 <- merge(pol2.df, actd.df, by = "tx_id")

##-------------------------------------------------------## For genes

pol1.df <- dplyr::select(res.pol1.filtered, baseMean, gene_id, log2FoldChange, padj, gene_biotype, diffexpressed)
pol2.df <- dplyr::select(res.pol2.filtered, baseMean, gene_id, log2FoldChange, padj, gene_biotype, diffexpressed)
actd.df <- dplyr::select(res.actd.filtered, baseMean, gene_id, log2FoldChange, padj, gene_biotype, diffexpressed)

## filter counts on DEG data
plot.1 <- merge(pol1.df, pol2.df, by = "gene_id")
plot.2 <- merge(pol1.df, actd.df, by = "gene_id")
plot.3 <- merge(pol2.df, actd.df, by = "gene_id")

##-------------------------------------------------------## Quarter data

## Function to plot LFC vs LFC for gene data. Plot = data to plot
plotLFC <- function(plot, x.lab, y.lab){
  
  require(viridis)
  require(ggpubr)
  require(ggplot2)
  require(ggside)
  
  plot$cond <- NA
  ## Split data up into quarter
  plot$cond[plot$log2FoldChange.x >0 & plot$log2FoldChange.y < 0] <- "Q4"
  plot$cond[plot$log2FoldChange.x >0 & plot$log2FoldChange.y > 0] <- "Q3"
  plot$cond[plot$log2FoldChange.x < 0 & plot$log2FoldChange.y > 0] <- "Q2"
  plot$cond[plot$log2FoldChange.x < 0 & plot$log2FoldChange.y < 0] <- "Q1"
  
  pos <- position_jitter(width = 0.2, seed = 1)  
  lm <- lm(log2FoldChange.y ~ log2FoldChange.x, plot) ## Get regression line
  coeff<-coefficients(lm)           
  intercept<-coeff[1] 
  slope<- coeff[2] 
  
  # Bin size control + color palette
  p <- ggplot(plot, aes(x=log2FoldChange.x, y=log2FoldChange.y) ) +
    geom_hex(bins = 100) +scale_y_continuous(limits = c(-10,10))+scale_x_continuous(limits = c(-10,10))+ 
    scale_fill_viridis( option = "B")+geom_vline(xintercept=0, col="darkred", lty = "dashed") +
    theme(plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
    geom_hline(yintercept=0, col="darkred", lty = "dashed") + theme_classic()+
    xlab(x.lab)+ ylab(y.lab) +
    geom_abline(intercept = 0, slope = 1, col = "darkred", linetype="dashed", size=0.5 )
  p + geom_xsidedensity(aes(y=after_stat(density), alpha=0.9, xfill = diffexpressed.x), position = "dodge") +
    geom_ysidedensity(aes(x=after_stat(density), alpha=0.9, yfill = diffexpressed.y), position = "dodge")+
    scale_xfill_manual(values =  c("#000004","#781c6d", "#ed6925")) +
    scale_yfill_manual(values =  c("#000004","#781c6d", "#ed6925")) + theme(ggside.panel.scale = .2)+
    scale_xsidey_continuous() + scale_ysidex_continuous() 
  
  print(p)
  
}

plotLFC(plot.1, x.lab = "Pol 1 LFC", y.lab = "Pol 2 LFC")
plotLFC(plot.2, x.lab = "Pol 1 LFC", y.lab = "ActD LFC")
plotLFC(plot.3, x.lab = "Pol 2 LFC", y.lab = "ActD LFC")

#########################################################################################
# SI Figure 4A-C: MA Plot using functions written for DESeq2
#########################################################################################

## Find coefficients for following function
resultsNames(dds)

### Call function for results
plot_MA("Conditions_ActD_6hrs_vs_CTRL_0hrs", "Act D", res.actd) ## ma.plot.actd
plot_MA("Conditions_Pol1_6hrs_vs_CTRL_0hrs", "Pol 1", res.pol1) ## ma.plot.pol1
plot_MA("Conditions_Pol2_6hrs_vs_CTRL_0hrs", "Pol 2", res.pol2)

#########################################################################################
# Base R Volacano plot of DESeq results
#########################################################################################

### Call function for results
DE_Vol_Plot(res.actd,"ActD: Differentially Expressed Genes", y.limits=300) ## volplot.actd.labels
DE_Vol_Plot(res.pol1,"Pol1: Differentially Expressed Genes", y.limits=300) ## volplot.pol1.labels
DE_Vol_Plot(res.pol2,"Pol2: Differentially Expressed Genes",y.limits=300) ## volplot.pol2.labels

### No labels (must comment out geom_text_repel() and label = ..)

#####################################################################################
# Bring in Loop, TAD, or NAD data here as genomic ranges                                     
#####################################################################################

library(plyranges)
library(GenomicRanges)
library(tidyverse) 

exps <- c("WT_HCT116_CTRL", "Pol1_6hrs_Aux", "Pol2_6hrs_Aux")

##-------------------------------------------------------## Loop BED
# (POLR1A Loops: N= 13,282; POLR2A Loops: N =56,420; DMSO Loops: N=18,123)
dir.loop <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/112123_HiC/112123_HiC/loop_analysis/"
list.files(dir.loop)
exps[3]
exps[1]
## Convert results to Genomic Ranges data
loops <- as.data.frame(read.table(file.path(dir.loop,exps[3], '/mega/hiccups_results/',"merged_loops.bedpe"), sep= "\t", header= F))
loops <- loops %>% 
  transform( seqnames= paste0("chr",loops$V1), start = loops$V22, end = loops$V23)  %>% 
  as_granges()

##-------------------------------------------------------## TAD BED

dir.tad <- '/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/112123_HiC/112123_HiC/contact_domains/arrowhead_domains'
list.files(dir.tad)
exps[2]

## Convert results to Genomic Ranges data
tads <- as.data.frame(read.table(file.path(dir.tad,exps[2],'/mega/inter_30_contact_domains', "5000_blocks.bedpe"), sep= "\t", header= F))
tads <- tads %>% 
  transform( seqnames= paste0("chr",tads$V1), start = tads$V2, end = tads$V3)  %>% 
  as_granges()

##-------------------------------------------------------## NAD BED

bed.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/NAD_tracks"
list.files(bed.dir)

nads <- as.data.frame(read_bed(file.path(bed.dir, "4DNFI4HQPGVC.bed") , col_names = NULL, genome_info = "hg38"))
nads

## Convert results to Genomic Ranges data
nads <- nads%>% 
  transform(seqnames, start, end)  %>% 
  as_granges()

#########################################################################################
# Density correlation between genes on same chromosome
#########################################################################################

## Prepare data
res.pol1.filtered <- res.pol1.filtered[res.pol1.filtered$diffexpressed == "Upregulated" | res.pol1.filtered$diffexpressed == "Downregulated",]
res.pol1.sorted<- res.pol1.filtered[order(res.pol1.filtered$seq_name, res.pol1.filtered$gene_seq_start),]
res.pol1.sorted <- mutate(res.pol1.sorted, gene.length = res.pol1.sorted$gene_seq_end - res.pol1.sorted$gene_seq_start)

res.pol2.filtered <- res.pol2.filtered[res.pol2.filtered$diffexpressed == "Upregulated" | res.pol2.filtered$diffexpressed == "Downregulated",]
res.pol2.sorted<- res.pol2.filtered[order(res.pol2.filtered$seq_name, -res.pol2.filtered$gene_seq_start),]
res.pol2.sorted <- mutate(res.pol2.sorted, gene.length = res.pol2.sorted$gene_seq_end - res.pol2.sorted$gene_seq_start)


## Count distance of nearest neighbors for each gene
res <- res.pol1.sorted$
res
tbl <- data.frame()
for (i in 1:nrow(res)){
  
  x <- res$gene_seq_start[i] + (res$gene.length[i]/2)
  y <- res$gene_seq_start[i+1] + (res$gene.length[i+1]/2)
  nn.dist = abs(y-x)
 
  
  df <- data.frame(res$seq_name[i], x, y, nn.dist, res$gene_biotype[i], res$log2FoldChange[i], res$gene_id[i])
  print(df)
  tbl <- rbind(tbl,df)
  
}

colnames(tbl) <- c("chrom", "start", "end", "nn.distance", "biotype", "lfc", "geneid")
tbl <- tbl[tbl$biotype == "protein_coding" | tbl$biotype == "lncRNA" ,]
tbl <- tbl[tbl$nn.distance < 5e06,]

#labels <- c("transcribed_processed_pseudogene" , "protein_coding", "TEC","transcribed_unitary_pseudogene", "processed_pseudogene", "unprocessed_pseudogene", "pseudogene")
#feat_tbl <- feat_tbl[feat_tbl$tx_type %in%labels ,]

## hist of values
colors <- c("#009392","#cf597e")
ggplot(tbl, aes(x=nn.distance, color = biotype, fill = "white")) + scale_color_brewer(palette="Dark2") + theme_classic()+
  geom_histogram(position="dodge", alpha = 1, bins=75, fill="white") +scale_x_continuous(limits=c(0,5e6))

## Density of values
tbl.1 <- mutate( tbl, name = paste0(tbl$biotype,"_Pol1"))
#tbl.2 <- mutate( tbl, name = paste0(tbl$biotype,"_Pol2"))

tbl <- rbind(tbl.1, tbl.2)
tbl <- rbind(tbl.1)

library(plyr)
mu <- ddply(tbl, "name", summarise, grp.mean=mean(nn.distance))
mu <- mu[-c(3:4),]
head(mu)

ggplot(tbl, aes(x=nn.distance, color = name, fill=name))  + theme_classic()+
  scale_color_brewer(palette="Dark2")+ scale_fill_brewer(palette="Dark2")+
  geom_density(alpha=0.2) + geom_vline(data=mu, aes(xintercept=grp.mean, color=name),
                              linetype="dashed") + theme(legend.position="bottom") 


## Plot each chromosome individually 
### make a null distribution and plot it along side (randomly grab 100 genes per chrom and calculate their nn distances)
tbl$factor<- as.numeric(ifelse(tbl$chrom == "X" , 23, tbl$chrom))

tbl <- data.frame(tbl)

plot <- tbl[tbl$name == "lncRNA_Pol1" | tbl$name == "protein_coding_Pol1",]
plot <- tbl[tbl$name == "lncRNA_Pol2" | tbl$name == "protein_coding_Pol2",]

plot %>%
  #mutate(chrom = fct_reorder(chrom, as.numeric(factor))) %>% 
  ggplot( aes(x=chrom, y=nn.distance, fill=name)) +
  geom_boxplot(notch=F, notchwidth = 0.8, outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_classic()
  #scale_y_continuous(labels = c("0 kb", "100 kb", "200 kb", "300 kb", "400 kb", "500 kb", "600 kb", "700 kb", "800 kb", "900 kb", "1 mb"),
                     #breaks = c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5, 8e5, 9e5, 1e6),
                     #limits = c(0, 1e6))+
  #theme(legend.position="none", plot.title = element_text(size=11), axis.text.x = element_text(angle = 90)) +
  #ggtitle("TADs remaining after Rad21 degradation") + xlab("") + ylab("Fraction of N")

#####################################################################################
# Pol 1 HiC Loop or TAD BED DEG intersection                                     
#####################################################################################

# note to self: 
# Try local S after that
# Try genes at anchors only after that (within 50 kb?)
# check CTCF and RAD21 tacks against loops and TADs
# check NADs for loops and TADs

##-------------------------------------------------------## Intro functions and packages

## Color scales
library(rcartocolor)
display_carto_all(type = 'diverging', colorblind_friendly = TRUE) ## get colorblind carto colors
my_colors = carto_pal(7, "Earth") 
my_colors = carto_pal(7, "ArmyRose") ## list carto colors
my_colors

## control: "#A16928" treatment: "#2887A1" ## Earth
## control: "#798234" treatment: "#D46780" ## ArmyRose

library(plyranges)
library(GenomicRanges)
library(tidyverse) 

write_out <- function(results, filename, dir= "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/results/"){
  write.csv(as.data.frame(results), file=paste0(dir,"/", filename, ".csv"))
}


## Declare experiments here
exps <- c("WT_HCT116_CTRL", "Pol1_6hrs_Aux")

dir.supp <- '/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/supp_data'
list.files(dir.supp)

## These are parsed loop and TAD files 
loops <- read.table(file.path(dir.supp, "pol1.loop.scaling.csv"), sep= ",", header= T) ## contact behavior for all loops
loops <- read.table(file.path(dir.supp, "pol1.tad.scaling.csv"), sep= ",", header= T) ## contact behavior for all tads
#loops <- read.table(file.path(dir.supp, "pol1.dge.scaling.inloops.csv"), sep= ",", header= T) ## contact behavior for genes inside loops
#loops <- read.table(file.path(dir.supp, "pol1.dge.scaling.intads.csv"), sep= ",", header= T) ##  contact behavior for genes inside tads

##-------------------------------------------------------##  Normalize for sequencing depth
## Use one or the other ## Do I need this? Matrix normalization should normalize data for sequence depth

norm.tbl <- loops %>% group_by(chrom,cond) %>% summarise( sum = mean(sum.count)) ## using the sum of observed counts
norm.tbl <- loops %>% group_by(chrom,cond) %>% summarise( sum = sum(mean.contact.prob)) ## Using the sum of mean contact probabilitie s

norm.tbl$chrom<- ifelse(norm.tbl$chrom == "X", as.numeric(23), as.numeric(norm.tbl$chrom))
loops$chrom<- ifelse(loops$chrom == "X", as.numeric(23), as.numeric(loops$chrom))

## Order table
norm.tbl <- norm.tbl[order(norm.tbl$chrom),]

## separate by condition
norm.tbl.wt <- norm.tbl[norm.tbl$cond== exps[1],]
norm.tbl.p1 <- norm.tbl[norm.tbl$cond== exps[2],]

## add in inverse condition's sums
norm.tbl.wt$sum.inv <- norm.tbl.p1$sum
norm.tbl.p1$sum.inv <- norm.tbl.wt$sum

## Bind and join norm table to loops table
norm.tbl <- rbind(norm.tbl.wt, norm.tbl.p1)
loops <- left_join(x = loops, y = norm.tbl, by = c("chrom", "cond"))

## Normalization step
loops <- loops %>% mutate(norm.counts=mean.counts/sum.inv) ## Normalize mean counts 
loops <- loops %>% mutate(norm.cont.prob=mean.contact.prob/sum.inv) ## Normalize mean contact behavior 
loops$chrom<- ifelse(loops$chrom == 23, "X", loops$chrom) ## replace chrom strings

##-------------------------------------------------------## merge back together

loops <- loops %>% 
  transform( seqnames= paste0("chr",loops$chrom), start = loops$start, end = loops$end)  %>% 
  as_granges()

##-------------------------------------------------------## For using counts instead of LFC
## If I just want the log normalized contacts, I use this

## grab normalized counts
pol1.counts <- rowMeans(assay(dds)[,4:6])
wt.counts <- rowMeans(assay(dds)[,10:12])
counts <- data.frame(cbind(pol1.counts, wt.counts), gene_id = names(wt.counts))

res.gr <- res.pol1.filtered[res.pol1.filtered$diffexpressed == "Upregulated" | res.pol1.filtered$diffexpressed == "Downregulated" ,] 

## filter counts on DEG data
counts <- counts[rownames(counts) %in% res.gr$gene_id,]
res.gr <- merge(counts, res.gr, by = "gene_id")

##-------------------------------------------------------## gene BED
## Grab deseq2 data for upregulated | downregulated genes

res.gr <- res.pol1.filtered[res.pol1.filtered$diffexpressed == "Upregulated" | res.pol1.filtered$diffexpressed == "Downregulated" ,] 
#res.gr <- res.gr[abs(res.gr$log2FoldChange) > 2,] ## Filter just high DEGs
res.gr<- res.gr %>% 
  transform( seqnames= paste0("chr",res.gr$seq_name), start = res.gr$gene_seq_start, end = res.gr$gene_seq_end)  %>% 
  as_granges()

##-------------------------------------------------------## NOTE
##
## Note: Use join_overlap_inner if you want all loop coordinates that 
## overlap with all gene coordinates. Use join_overlap_intersect if you
## want gene coordinates that overlap with loop coordinates i.e. loop
## chunks. join_overlap_intersect will merge adjacent or overlapping genes
## FYI. Both are valid. Use x=loops, y=genes
##
##-------------------------------------------------------## filter non Overlaps

## get list of all loops that overlap with genes
overlaps <- data.frame(join_overlap_inner(loops,res.gr)) 

##-------------------------------------------------------## Write out gene overlaps here
## this section is just for getting the genes that overlap with their full coordinates

## prepare gene data
exp <- exps[2]
gr <- loops %>% filter(loops$cond==exp)

## Count overlaps
tbl.out <- data.frame(res.gr %>%
        mutate(n_overlaps = count_overlaps(., gr, maxgap = 0)) %>% 
        filter(n_overlaps > 1) %>% mutate(name = rep(exp, length(n_overlaps > 1))))

## Get distinct genes for each conditions
tbl.wt <- distinct(tbl.out, tbl.out$gene_id, .keep_all = T)
tbl.p1 <- distinct(tbl.out, tbl.out$gene_id, .keep_all = T)
tbl.out <- rbind(tbl.wt, tbl.p1)

## Write out
write_out(tbl.out, "all.dge.loops.plyranges")
write_out(tbl.out, "all.dge.tads.plyranges")

##-------------------------------------------------------## Plot merged data
## Here, plot the merged counts and DEG data

## control: "#A16928" treatment: "#2887A1" ## Earth
## control: "#798234" treatment: "#D46780" ## ArmyRose

library(ggplot2)
library("ggExtra")
library("ggside")
library(PupillometryR)

## Separate data and get distinct values
tbl.1 <- overlaps[overlaps$cond == exps[1],] %>%  distinct(gene_id, .keep_all = T)
tbl.2 <- overlaps[overlaps$cond == exps[2],] %>%  distinct(gene_id, .keep_all = T)
tbl <- rbind(tbl.1, tbl.2) ## merge back together

# filename: loops.lfc.vs.cp.dotplot | tads.lfc.vs.cp.dotplot
p <- ggplot(overlaps, aes(x=log2FoldChange, y=log10(mean.contact.prob), color=cond)) +theme_classic() +
  scale_color_manual(values=c("#D46780","#798234"))+
  scale_x_continuous(limits = c(-5,5))+ 
  scale_y_continuous(limits = c(-3.5,-1))+ 
  geom_point(size=1,alpha = 0.6 ) + theme(legend.position="bottom")+
  theme(legend.position="none", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("WT & Pol1 TADs | Contact Probability vs DEGs") + xlab("Log Fold Change") + ylab("Mean Contact Probability")
p

# filename: loops.mean.cp.rainplot | tad.mean.cp.rainplot | loops.mean.cnts.rainplot | tad.mean.cnts.rainplot
overlaps %>%
  #ggplot( aes(x=cond, y=log10(mean.contact.prob), fill=cond)) +
  ggplot( aes(x=cond, y=mean.counts, fill=cond)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15))+ coord_flip() +
  scale_fill_manual(values=c("#D46780","#798234"))+
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
  theme_classic() + geom_boxplot(width = 0.2,notch=TRUE, notchwidth = 0.3, outlier.shape = NA, color = "white")+
  #scale_y_continuous(labels = c("-0.5","-0.4", "-0.3","-0.2", "-0.1", "0"), breaks = c(-0.5,-0.4, -0.3,-0.2, -0.1, 0), limits = c(-0.5,0))+ 
  #scale_y_continuous(limits = c(-3.5,-1))+ 
  theme(legend.position="none", plot.title = element_text(size=14), text = element_text(size=14, family="Arial"), axis.text.x = element_text(angle = 0)) +
  #ggtitle("WT & Pol1 TADs | Contact Probability ") + xlab("") + ylab("Mean Contact Probability")
  ggtitle("WT & Pol1 TADs | Observed Counts ") + xlab("") + ylab("Log10 Mean Observed Counts")

##-------------------------------------------------------## Plot count data

## Generate regions data to filter on 
overlaps$regions = paste0(overlaps$seqnames,'.',overlaps$start,'.',overlaps$end, ".", overlaps$cond) ## get regions for filtering
length(unique(overlaps$regions)) 

tbl <-overlaps%>%
  reframe(down.lfc = sum(down.lfc), up.lfc = sum(up.lfc),down = sum(down), up = sum(up),
          genes =list(unique(symbol)), biotype =list(unique(gene_biotype)), slope = unique(unique(slope)), 
          Pol1_6hrs_Aux = mean(pol1.counts), WT_HCT116_CTRL= mean(wt.counts), ## for count data
          mean.contact.prob=unique(mean.contact.prob), mean.counts = unique(mean.counts),.by = regions) %>% 
  separate_wider_delim(regions, ".", names = c("chr", "start", "end", "cond"))

## Gather count data
counts <- gather(tbl, key="name",value = "counts", 12:13,)

## Filter
counts.wt <- counts %>% filter(cond == exps[1] & name == exps[1])
counts.p1 <- counts %>% filter(cond == exps[2] & name == exps[2])
counts <- rbind(counts.wt, counts.p1)

## filename: loops.cp.vs.rawcts.margplot | tads.cp.vs.rawcts.margplot
p <- ggplot(counts, aes(x=log10(mean.contact.prob), y=log10(counts), color=cond)) +theme_classic() +
  scale_color_manual(values=c("#D46780","#798234"))+
  geom_point(size=1,alpha = 0.8 ) + theme(legend.position="bottom") + #scale_y_continuous(limits = c(1,1e7))+
  theme( plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("WT & Pol1 TADs | Raw counts vs CP ") + xlab("Raw Counts") + ylab("Mean Contact Probability ")
p + geom_xsidedensity(aes(y=after_stat(density), alpha=0.5, xfill = cond), position = "dodge") +
  geom_ysidedensity(aes(x=after_stat(density), alpha=0.5, yfill = cond), position = "dodge")+
  scale_xfill_manual(values = c("#D46780","#798234")) +
  scale_yfill_manual(values = c("#D46780","#798234")) + theme(ggside.panel.scale = .2)+
  scale_xsidey_continuous() + scale_ysidex_continuous() 

counts  %>% # loops.rawcts.rainplot | tads.rawcts.rainplot 
  ggplot( aes(x=name, y=log10(counts), fill=name)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15))+ coord_flip() +
  scale_fill_manual(values=c("#D46780","#798234"))+
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
  theme_classic() + geom_boxplot(width = 0.2,notch=TRUE, notchwidth = 0.3, outlier.shape = NA, color = "white")+
  #scale_y_continuous(labels = c("-0.5","-0.4", "-0.3","-0.2", "-0.1", "0"), breaks = c(-0.5,-0.4, -0.3,-0.2, -0.1, 0), limits = c(-0.5,0))+ 
  #scale_y_continuous(limits = c(0,0.02))+ 
  theme(legend.position="none", plot.title = element_text(size=11), axis.text.x = element_text(angle = 0)) +
  ggtitle("Raw counts in TADs") + xlab("")

##-------------------------------------------------------## Plot CDFs

## Gather and count number of log2 points 
ko_tbl <- overlaps[overlaps$cond ==exps[2],]
ctrl_tbl <- overlaps[overlaps$cond ==exps[1],]

## CDF for ko vs ctrl
ko.cdf <- CDF_X(ko_tbl$mean.contact.prob, "ko")
ctrl.cdf <- CDF_X(ctrl_tbl$mean.contact.prob, "ctrl")

## Bind
cdf.plot <- rbind(ko.cdf,ctrl.cdf)

## Plot CDF ## filename: loops.cp.cdfplot | tads.cp.cdfplot
ggplot(cdf.plot, aes(x=values, y=cumsums, group = name, color = name)) +
  geom_line() + theme_classic()+
  #scale_y_continuous(limits = c(0,  1), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
  scale_x_continuous(limits = c(0,0.03))+
  scale_color_manual(values=c("#798234","#D46780")) + 
  theme( plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("WT & Pol1 tads | Contact Probability ") + xlab("CDF") + ylab("Mean Contact Probability")

#####################################################################################
# NAD BED DEG | Loop | TAD intersection                                     
#####################################################################################

##-------------------------------------------------------## NAD overlap

## Find intersection between LADs and all genes
nad.overlaps <- data.frame(join_overlap_inner(nads, res.gr))

## Find intersection between LADs and all loops or tads
nad.overlaps <- data.frame(join_overlap_inner(nads, loops))

## Find non-overlapping regions of genes
non.nads<- data.frame(res.gr[!res.gr %over% nad.overlaps,])

## Find non-overlapping regions of loops
non.nads<- data.frame(loops[!loops %over% nads,])

##----------------------------------NAD DEG data exploration----------------------------------##

## get table of number of DEGs up and down in  and outside of NADs
feat_tbl.1 <- data.frame(nad.overlaps) %>% group_by(gene_biotype) %>% 
  summarise(down = -sum(down), up = sum(up))
feat_tbl.2 <- data.frame(non.nads) %>% group_by(gene_biotype) %>% 
  summarise(down = -sum(down), up = sum(up))

## get table of magnitude of DEGs up and down (LFC) in and outside of NADs
feat_tbl.1 <- data.frame(nad.overlaps) %>% group_by(gene_biotype) %>% 
  summarise(down = sum(down.lfc), up = sum(up.lfc))
feat_tbl.2 <- res.data.frame(non.nads) %>% group_by(gene_biotype) %>% 
  summarise(down = sum(down.lfc), up = sum(up.lfc))

##----------------------------------NAD DEG data exploration----------------------------------##
library(ggridges)

nad.overlaps <- nad.overlaps %>% mutate(cond = paste0(nad.overlaps$cond,"_", "nad_ol"))
non.nads <- non.nads %>% mutate(cond = paste0(non.nads$cond,"_", "nad_non_ol"))

plot.nads <- rbind(non.nads, nad.overlaps)

## Color scales
colors <- colorRampPalette(c("white","#798234","#D46780"))(8)
colors <- c("#D8DBC5", "#A5AB77", "#97794D", "#BA6E6A")

## Plot NAD CP ## filename: loops.nads.cp.distplot | tads.nads.cp.distplot | loops.nads.cp.log10.distplot | tads.nads.cp.log10.distplot
ggplot(plot.nads, aes(x = log10(mean.contact.prob), y = cond, fill = stat(quantile))) +
  stat_density_ridges(quantile_lines = F,calc_ecdf = TRUE,geom = "density_ridges_gradient") +theme_classic() +
  scale_y_discrete(breaks = c("WT_HCT116_CTRL_nad_ol", "WT_HCT116_CTRL_nad_non_ol", "Pol1_6hrs_Aux_nad_ol", "Pol1_6hrs_Aux_nad_non_ol"),
                   labels = c("WT NAD OL", "WT NAD Non-OL", "Pol1 NAD OL", "Pol1 NAD Non-OL"))+
  scale_fill_manual(values=colors) + theme( legend.position="none",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("TAD Contact Probability in/outside of NADs | Loops") + xlab("Log 10 Mean Contact Probability") + ylab("")


#########################################################################################
# Plotting TAD scaling vs TAD size For Pol 1
#########################################################################################
library(ggplot2)
library("ggExtra")
library(ggplot2)
library(PupillometryR)
library("ggside")

dir.supp <- '/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/supp_data'
list.files(dir.supp)

tad.s <- read.table(file.path(dir.supp, "pol1.tad.scaling.csv"), sep= ",", header= T) ## tads
tad.s <- read.table(file.path(dir.supp, "pol1.loop.scaling.csv"), sep= ",", header= T) ## loops

tad.s <- mutate(tad.s, size = end-start)

## filename: loops.size.vs.s.margplot| tads.size.vs.s.margplot
p <- ggplot(tad.s, aes(x=tad.s$slope, y=tad.s$log2FoldChange, color=cond)) +theme_classic() +
  scale_color_manual(values=c("#D46780","#798234"))+
  scale_y_continuous(labels = c("0", "0.25 mbp", "0.50 mbp", "0.75 mbp", "1.00 mbp", "1.25 mbp", "1.50 mbp", "1.75 mbp", "2.00 mbp","2.25 mbp", "2.50 mbp", "2.75 mbp", "3.00 mbp"),breaks = c(0, 2.5e5, 5e5, 7.5e5, 1e6, 1.25e6, 1.5e6, 1.75e6, 2e6, 2.25e6, 2.5e6, 2.75e6, 3e6),limits = c(0,3e6))+ 
  scale_x_continuous(labels = c("-0.6", "-0.5","-0.4", "-0.3","-0.2", "-0.1", "0"), breaks = c(-0.6, -0.5,-0.4, -0.3,-0.2, -0.1, 0),limits = c(-0.5,0))+ 
  geom_point(size=1,alpha = 0.8 ) + theme(legend.position="bottom") + #scale_y_continuous(limits = c(1,1e7))+
  theme( plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("WT & Pol1 TADs | Scaling vs Size ") + xlab("Contact Scaling") + ylab("Size")
p + geom_xsidedensity(aes(y=after_stat(density), alpha=0.5, xfill = cond), position = "dodge") +
  geom_ysidedensity(aes(x=after_stat(density), alpha=0.5, yfill = cond), position = "dodge")+
  scale_xfill_manual(values = c("#D46780","#798234")) +
  scale_yfill_manual(values = c("#D46780","#798234")) + theme(ggside.panel.scale = .2)+
  scale_xsidey_continuous() + scale_ysidex_continuous() 



#########################################################################################
# Hg38 gene positions
#########################################################################################

#chrom.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Raw_Data/Annotations"
chrom.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations"
list.files(chrom.dir)

## Genomic coordinates for all coding genes in hg38
gene.pos <- read.table(file.path(chrom.dir, "hg38_gene_positions.csv"), sep= ",", header= T)

## Convert to Genomic Ranges data
gr.gene.pos <- as.data.frame(gene.pos)
gr.gene.pos <- gr.gene.pos %>% 
  transform( seqnames= paste0("chr",gr.gene.pos$chromosome_name), start = gr.gene.pos$start_position, end = gr.gene.pos$end_position)  %>% 
  as_granges()
gr.gene.pos


#########################################################################################
# Bin NAD and Non-NAD areas of the genome
#########################################################################################

library(plyranges)
library(dplyr)
library(tidyr)

## Use this for correlation analysis too ##
## Excellent guide on correlation analysis in R: http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r

exps <- c("WT_HCT116_CTRL","Pol1_6hrs_Aux")

## gene position directory
chrom.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations"
list.files(chrom.dir)

## chrom sizes for mm10
chrom.sizes <- read.table(file.path(chrom.dir, "hg38.chrom.sizes"), sep= "\t", header= F)
colnames(chrom.sizes) <- c("name", "start", "size")

## For loop to generate 1 MB segments of the genome for counting overlaps in peaks
tbl <- data.frame()
bin.size = 2000000 #Bins = 10k, 25k, 100k, and 1 MB
for (i in 1:length(chrom.sizes$name)){
  
  chrom <- chrom.sizes$name[i]
  
  x <-as.numeric(seq(0,chrom.sizes$size[i], by = bin.size))
  
  y <-as.numeric(seq(bin.size,chrom.sizes$size[i], by = bin.size))
  y <- append(y, chrom.sizes$size[i], after = length(y))
  
  chr <- rep(chrom, length(x))
  cond <- rep(exps[1], length(x))
  size = y-x
  
  df <- data.frame(chr, start=x, end=y, name=cond, size = size)
  print(df)
  tbl <- rbind(tbl,df)
  
}

hg38.tbl.wt <- tbl
hg38.tbl.p1 <- tbl

hg38.tbl <- rbind(hg38.tbl.wt , hg38.tbl.p1)

## Generate ranges from above
hg38.gr <- hg38.tbl %>% 
  transform( seqnames= chr, start = start, end = end)  %>% 
  as_granges()

##-------------------------------------------------------## use NAD BED instead

## Bin NADs for S calculation
tbl <- data.frame()
bin.size = 100000 #Bins = 10k, 25k, 100k, and 1 MB # 5 mb doesn't work (too big)
for (i in 1:length(nads$seqnames)){
  
  chrom <- nads$seqnames[i]
  
  x <-as.numeric(seq(nads$start[i]-1,nads$end[i]-bin.size, by = bin.size))
  
  y <-as.numeric(seq(nads$start[i]-1 +bin.size,nads$end[i], by = bin.size))
  #y <- append(y, nads$width[i], after = length(y))
  
  chr <- rep(chrom, length(x))
  cond <- rep(exps[1], length(x))
  size = y-x
  
  df <- data.frame(chr, start=x, end=y, name=cond, size = size)
  print(df)
  
  tbl <- rbind(tbl,df)
  
}

##-------------------------------------------------------## process whole NADs

nad.tbl.wt <- nads %>% mutate(strand = rep(exps[1], nrow(nads)))%>% rename(size = width, name = strand, chr = seqnames) 
nad.tbl.p1 <- nads %>% mutate(strand = rep(exps[2], nrow(nads)))%>% rename(size = width, name = strand, chr = seqnames) 

##-------------------------------------------------------## process whole NADs
nad.tbl.wt <- tbl
nad.tbl.p1 <- tbl

nad.tbl <- rbind(nad.tbl.p1, nad.tbl.wt)

## Convert results to Genomic Ranges data
nads.gr <- nad.tbl %>% 
  transform( seqnames= chr, start = start, end = end)  %>% 
  as_granges()

##-------------------------------------------------------## NAD GR overlap

## Find non-overlapping regions of genes
non.nads <- data.frame(hg38.gr[!hg38.gr %over% nads.gr,])

##-------------------------------------------------------## Prep data and write out
dir.out <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations"

## Clean up data
nad.tbl <- nad.tbl %>% mutate(cond = paste0("nad")) %>% mutate(chr=sapply(strsplit(as.character(chr), "chr"), `[`, 2))
non.nads <- non.nads %>% mutate(cond = paste0("non_nad")) %>% mutate(seqnames=sapply(strsplit(as.character(seqnames), "chr"), `[`, 2)) %>% 
                         select(-width,-strand, -chr) %>% rename(chr = seqnames)

## Write out
nad.tbl.out <- rbind(non.nads,nad.tbl)
write.csv(nad.tbl.out , file=paste0(dir.out,"/whole.nad.positions.csv"), row.names=FALSE)


##-------------------------------------------------------## Plot binned NAD Scaling Results
library(ggridges)
library(ggplot2)
library(forcats)

res.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/results/local_scaling_results"
list.files(res.dir)

## chrom sizes for mm10
scaling<- read.table(file.path(res.dir, "whole.nad.scaling.csv"), sep= ",", header= T)
colnames(scaling) <- c("cond", "chrom", "pos", "start", "end","log10start", "log10end", "slope", "intercept", "mean.contact.prob", "mean.counts")
plot.nads<- scaling%>% mutate(cond = paste0(cond,"_", pos))

## Color scales
colors <- colorRampPalette(c("white","#798234","#D46780"))(8)
colors <- c("#D8DBC5", "#A5AB77", "#97794D", "#BA6E6A")


## Plot NAD CP ## filename: whole.nad.scaling.1e5to1e6.dist | whole.nad.cp.1e5to1e6.dist | whole.nad.scaling.dist | whole.nad.cp.dist 
plot.nads %>% mutate(cond = fct_relevel(cond, c("Pol1_6hrs_Aux_non_nad","WT_HCT116_CTRL_non_nad", "Pol1_6hrs_Aux_nad", "WT_HCT116_CTRL_nad" ))) %>% 
ggplot(aes(x = slope, y = cond, fill = stat(quantile))) +
  stat_density_ridges(quantile_lines = F,calc_ecdf = TRUE,geom = "density_ridges_gradient") +theme_classic() +
  scale_y_discrete(breaks = c("Pol1_6hrs_Aux_non_nad","WT_HCT116_CTRL_non_nad", "Pol1_6hrs_Aux_nad", "WT_HCT116_CTRL_nad" ),
                   labels = c("Pol1 non-NADs", "WT non-NADs","Pol1 NADs" , "WT NADs"))+
  #scale_x_continuous(limits = c(0, 0.00005))+
  scale_fill_manual(values=colors) + theme( legend.position="none",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("Contact Scaling in/outside of NADs | All Length Scales") + xlab("Contact Scaling") + ylab("")
  #xlab("log10 Mean Contact Probability") + ylab("")


#########################################################################################
# Nearest Neighbor Analysis of DEGs vs ChIP marks
#########################################################################################

bed.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations"
list.files(bed.dir)

h3k9.bed <- read.table(file.path(bed.dir,"/H3K9Me3/", "ENCFF832IOO.bed"), sep= "\t", header= F)
h3k4.bed <- read.table(file.path(bed.dir,"/H3K4me3/", "ENCFF711MPL.bed"), sep= "\t", header= F)
ctcf.bed <- read.table(file.path(bed.dir,"/ctcf/", "ENCFF463FGL.bed"), sep= "\t", header= F)
rad21.bed <- read.table(file.path(bed.dir,"/rad21/", "ENCFF391AAM.bed"), sep= "\t", header= F)
pol2.bed <- read.table(file.path(bed.dir,"/polr2a/", "ENCFF848IHI.bed"), sep= "\t", header= F)
pol2ps5.bed <- read.table(file.path(bed.dir,"/polr2a-ps5/", "ENCFF229YMU.bed"), sep= "\t", header= F)


nndist <- function(res, bed, diffexp){
  
require(plyranges)
  
  if ( diffexp == TRUE) {
    
    res <- res[res$diffexpressed == "Upregulated" | res$diffexpressed == "Downregulated" ,] 
    #res <- res[res$diffexpressed == "Downregulated" ,] 
    #res <- res[res$diffexpressed == "Upregulated" ,] 
    res.gr<- res %>% 
      #transform( seqnames= paste0("chr",res$seq_name), start = res$gene_seq_start, end = res$gene_seq_end)  %>% 
      transform( seqnames= paste0("chr",res$seq_name), start = res$tx_seq_start, end = res$tx_seq_end)  %>% 
      as_granges()
  } else {
    
    res <- res[res$diffexpressed != "Upregulated" & res$diffexpressed != "Downregulated" ,] 
    #res <- res[res$diffexpressed != "Downregulated",] 
    #res <- res[res$diffexpressed != "Upregulated",] 
    res.gr<- res %>% 
      #transform( seqnames= paste0("chr",res$seq_name), start = res$gene_seq_start, end = res$gene_seq_end)  %>% 
      transform( seqnames= paste0("chr",res$seq_name), start = res$tx_seq_start, end = res$tx_seq_end)  %>% 
      as_granges()
  }

    bed.f <- bed %>% 
      transform( seqnames= paste0(bed$V1), start = bed$V2, end = bed$V3)  %>% 
      as_granges()

    x <- data.frame(plyranges::add_nearest_distance(res.gr , bed.f, name = "distance"))
    return(x)

}

## h3k9
dge.pol1.nn <- nndist(res.pol1.filtered, h3k9.bed , TRUE)
pol1.nn <- nndist(res.pol1.filtered, h3k9.bed , FALSE)

## h3k4
dge.pol1.nn <- nndist(res.pol1.filtered, h3k4.bed , TRUE)
pol1.nn <- nndist(res.pol1.filtered, h3k4.bed , FALSE)

## ctcf
dge.pol1.nn <- nndist(res.pol1.filtered, ctcf.bed , TRUE)
pol1.nn <- nndist(res.pol1.filtered, ctcf.bed , FALSE)

## rad21
dge.pol1.nn <- nndist(res.pol1.filtered, rad21.bed , TRUE)
pol1.nn <- nndist(res.pol1.filtered, rad21.bed , FALSE)

## pol2
dge.pol1.nn <- nndist(res.pol1.filtered, pol2.bed , TRUE)
pol1.nn <- nndist(res.pol1.filtered, pol2.bed , FALSE)

## pol2-ps5
dge.pol1.nn <- nndist(res.pol1.filtered, pol2ps5.bed , TRUE)
pol1.nn <- nndist(res.pol1.filtered, pol2ps5.bed , FALSE)


##-------------------------------------------------------## Decile plots

## Found counts not that useful. Used BaseMean for all samples (baseMean uses all samples)
## And LFC for individual samples when looking at deciles

head(counts(dds, normalized=TRUE), n=1)

## Function for grabbing the baseMean normalized counts for specficic conditions
getSampleBM <- function(res, vec){

    cat( "Using counts from", colnames(head(counts(dds, normalized=TRUE)[,vec], n=1)))
    
    samp.baseMean <- rowMeans(counts(dds, normalized=TRUE)[,vec])
    
    res <- res[res$gene_id %in% names(samp.baseMean),]
    
    res$samp.baseMean <- samp.baseMean[names(samp.baseMean) %in% res$gene_id]
    
    return(res)

}

## Get counts
res.actd.filtered <- getSampleBM(res.actd.filtered, vec=1:3) ## ActD cond
res.pol1.filtered <- getSampleBM(res.pol1.filtered, vec=4:6) ## Pol 1 cond
res.pol1.filtered.wt <- getSampleBM(res.pol1.filtered, vec=10:12) ## WT cond

nndist_Tile <- function(res, bed, quantile, cond.label, diffexp = TRUE){
  
  require(plyranges)
  require(dplyr)
    
    if ( diffexp == TRUE) {
      
      #res <- res[res$diffexpressed == "Upregulated" | res$diffexpressed == "Downregulated" ,] 
      #res <- res[res$diffexpressed == "Downregulated" ,] 
      res <- res[res$diffexpressed == "Upregulated" ,] 
      res.gr<- res %>% 
        transform( seqnames= paste0("chr",res$seq_name), start = res$gene_seq_start, end = res$gene_seq_end)  %>% 
        as_granges()
    } else {
      
      res.gr<- res %>% 
        transform( seqnames= paste0("chr",res$seq_name), start = res$gene_seq_start, end = res$gene_seq_end)  %>% 
        as_granges()
    }

    bed.f <- bed %>% 
      transform( seqnames= paste0(bed$V1), start = bed$V2, end = bed$V3)  %>% 
      as_granges()
    
    
    x <- data.frame(plyranges::add_nearest_distance(res.gr , bed.f, name = "distance"))
    
    #x$bin <- ntile(x$baseMean, as.numeric(quantile))
    #x$bin <- ntile(x$samp.baseMean, as.numeric(quantile)) ## sampleBaseMeans instead of all counts
    x$bin <- ntile(abs(x$log2FoldChange ), as.numeric(quantile)) ## use for LFC instead of counts
    
    df <- x %>% group_by(bin) %>% reframe(mean.dist = log10(mean(na.omit(distance))), 
                                          mean.cnt = log10(mean(baseMean)), ## ALL Counts
                                          #mean.cnt = log10(mean(samp.baseMean)), ## Sample Counts
                                          mean.lfc = mean(abs(log2FoldChange)))
    df$cond <- rep(cond.label, nrow(df))
  
    return(df)
  
}

## Pick condition 
res = res.pol1.filtered
res = res.actd.filtered
res = res.pol1.filtered.wt

## Call function 
nn.h3k4 <- nndist_Tile(res, h3k4.bed , 10, "h3k4")
nn.h3k9 <- nndist_Tile(res, h3k9.bed , 10, "h3k9")

nn.ctcf <- nndist_Tile(res,ctcf.bed , 10, "CTCF")
nn.rad21 <- nndist_Tile(res, rad21.bed , 10, "Rad21")

nn.pol2 <- nndist_Tile(res,pol2.bed , 10, "Pol2")
nn.pol2ps5 <- nndist_Tile(res, pol2ps5.bed , 10, "Pol2-PS5")

## Bind and test plot 
df <- rbind(nn.h3k4, nn.h3k9)
df <- rbind(nn.rad21, nn.ctcf)
df <- rbind(nn.pol2, nn.pol2ps5)

## Test plot 
ggplot(df, aes(x=mean.lfc, y=mean.dist, color = cond)) + 
  geom_point()

colors <- c("#D8DBC5", "#A5AB77", "#97794D", "#BA6E6A")

## Full plot 
mark1 = "POL2"
mark2 = "POL2-PS5"

## h3k9.h3k4.decile.lineplot ## rad21.ctcf.decile.lineplot ## pol2.pol2ps5.decile.lineplot
ggplot(df, aes(x=mean.cnt, y=mean.dist, color = cond, fill = cond)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor( label.y.npc="bottom", label.x.npc = "left", r.digits = 2, method = "spearman") +
  scale_color_manual(values = c("#A5AB77", "#BA6E6A"), aesthetics = c("color", "fill")) +
  theme_pubr() + ggtitle(paste0("Gene Expression Vs Distance: ", mark1," & ", mark2))+
  labs(x = "Log10 Mean Expression", y = "Log10 Mean Distance")

##-------------------------------------------------------## Pheatmap of decile data

library(pheatmap)

## generate matrix
df.dist <- data.frame(H3K4 = nn.h3k4$mean.dist, H3K9 = nn.h3k9$mean.dist,RAD21 = nn.rad21$mean.dist, CTCF = nn.ctcf$mean.dist,POL2 = nn.pol2$mean.dist, POL2PS5 = nn.pol2ps5$mean.dist)
df.cnt <- nn.h3k4$mean.cnt ## get count vector
df.lfc <- nn.h3k4$mean.lfc ## get lfc vector

## Color scales
colors <- colorRampPalette(c("#798234","white","#D46780"))(64)

## Final Value – Initial Value) ÷ | Initial Value | × 100 = PERCENTAGE CHANGE.
rownames <- paste0(seq(1,10,1),": ", round(df.lfc, 2))

## Get heatmap of distance data ## decile.heatmap.allcounts ## decile.heatmap.pol1.lfc.up ## decile.heatmap.pol1.lfc.dn ## decile.heatmap.actd.lfc.dn ## decile.heatmap.actd.lfc.up
pheatmap(df.dist, display_numbers = TRUE, number_color = "black", color = colors,
         labels_row = rownames,
         fontsize_number = 8,  cluster_rows = FALSE,  cluster_cols = T, 
         #main = "Mean Distance Binned By Expression (All Counts)"
         main = "Mean Distance Binned By Act D Upregulated LFC"
         )

## use one column as vector for comparing mean log 10 counts
pheatmap(df.cnt, display_numbers = TRUE, number_color = "black", color = colors,
         fontsize_number = 8,  cluster_rows = F,  cluster_cols = F)
##-------------------------------------------------------## Plot CDFs

## CDF for ko vs ctrl
deg.cdf <- CDF_X(dge.pol1.nn$distance, "dge")
ctrl.cdf <- CDF_X(pol1.nn$distance, "all")

## Bind
cdf.plot <- rbind(deg.cdf,ctrl.cdf)

## Raw plot using gg
## Plot CDF ## filename: loops.cp.cdfplot | tads.cp.cdfplot
ggplot(cdf.plot, aes(x=values, y=cumsums, group = name, color = name)) +
  geom_line() + theme_classic()+
  #scale_y_continuous(limits = c(0,  1), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
  #scale_x_continuous(limits = c(0,1.5e6))+ ## K9
  scale_x_continuous(limits = c(0,1.5e6))+
  scale_color_manual(values=c("#798234","#D46780")) + 
  theme( plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("pol2 proximity to nearest gene") + xlab("distance in bp") + ylab("Mean Contact Probability")+   
  geom_segment(aes(x = mean(values[name=="dge"]), xend = mean(values[name=="dge"]), y = 0, yend = 1), linetype = "dashed")+
  geom_segment(aes(x = mean(values[name=="all"]), xend = mean(values[name=="all"]), y = 0, yend = 1), linetype = "dashed")

## Plot using ggpubr function (looks cleaner)
library(ggpubr)
library(plyr)

## Get stats for plotting mean
cdf.plot <- cdf.plot[cdf.plot$values < 1e6,]
cdat <- ddply(cdf.plot, "name", summarise, mean.v=mean(values))

label = "Pol 2 ps5"
limit.y = 1.5e6

## Use this function to plot CDF values ## h3k9.pol1.nncdf ## h3k4.pol1.nncdf ## ctcf.pol1.nncdf ## rad21.pol1.nncdf ## pol2.pol1.nncdf
ggecdf(cdf.plot, x="values", color = "name", linetype = "name",
       palette = c("#D46780", "#798234"), size = 1, ggtheme = theme_classic()) +
  scale_x_continuous(limits = c(0,limit.y ))+  
  geom_vline(data=cdat, aes(xintercept=mean.v,  colour=name),linetype="dashed", size=1)+ 
  theme( plot.title = element_text(size=12), text = element_text(size=16, family="Arial")) +
  labs(title=paste0("DEG Distance To Nearest ", label), 
       subtitle=paste0(label, ": Mean DEG Dist: ", as.integer(mean(cdf.plot$values[cdf.plot$name=="dge"])), 
                       " BP | Mean Gene Dist: ",as.integer(mean(cdf.plot$values[cdf.plot$name=="all"])), " BP"))+
                      xlab("Distance (BP)") + ylab("CDF")+
  scale_y_continuous(limits = c(0,  1), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  
##-------- ## set this code aside for later

# minimum value of x for the area under the curve shading
x_start <- deg.cdf[deg.cdf$cumsums ==.25,]
x_end   <- deg.cdf$cumsums <.50

#Subset the data and add the coordinates to make it shade to y = 0
shade <- rbind(c(x_start,0), subset(deg.cdf, values >= 
                                      x_start & values <= x_end), c(x_end, 0))

p <- ggplot(deg.cdf, aes(x=values, y=cumsums)) + 
  geom_line() + 
  xlab("Time") + 
  ylab("Cumulative Density")  + theme_classic()
p
shade
# add shading to cdf curve
p + geom_polygon(data = shade, aes(x=values , y=cumsums))

##-------- ## set this code aside for later

##-------------------------------------------------------## Plot dists

dge.pol1.nn$cond <- rep("deg", nrow(dge.pol1.nn))
pol1.nn$cond <- rep("all", nrow(pol1.nn))

data <- rbind(dge.pol1.nn, pol1.nn)
data$distance <- ifelse(data$distance ==0,log10(1), log10(data$distance))
data$baseMean <- log10(data$baseMean)

ggdensity(data, x = "distance",
          add = "mean", rug = TRUE,
          color = "cond", fill = "cond",
          palette = c("#00AFBB", "#E7B800"))


cols <- c("#F76D5E", "#FFFFBF", "#72D8FF")
# Density areas without lines
ggplot(data, aes(x = distance, fill = cond)) +
  geom_density(alpha = 0.8, color = NA) + 
  scale_fill_manual(values = cols) + theme_classic()



## Color scales
colors <- colorRampPalette(c("white","#798234","#D46780"))(8)
colors <- c("#D8DBC5", "#A5AB77", "#97794D", "#BA6E6A")


## Plot NAD CP ## filename: whole.nad.scaling.1e5to1e6.dist | whole.nad.cp.1e5to1e6.dist | whole.nad.scaling.dist | whole.nad.cp.dist 
data %>%
  ggplot(aes(x = distance, y = cond, fill = stat(quantile))) +
  stat_density_ridges(quantile_lines = F,calc_ecdf = TRUE,geom = "density_ridges_gradient") +theme_classic() +
  #scale_y_discrete(breaks = c("Pol1_6hrs_Aux_non_nad","WT_HCT116_CTRL_non_nad", "Pol1_6hrs_Aux_nad", "WT_HCT116_CTRL_nad" ),
                   #labels = c("Pol1 non-NADs", "WT non-NADs","Pol1 NADs" , "WT NADs"))+
  #scale_x_continuous(limits = c(0, 1e5))+
  scale_fill_manual(values=colors) + theme( legend.position="none",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("Contact Scaling in/outside of NADs | All Length Scales") + xlab("Contact Scaling") + ylab("")
#xlab("log10 Mean Contact Probability") + ylab("")

library(ggplot2)
data <- pol1.anno1
dat <- data[data$annotation == "Distal Intergenic",]

dat$distanceToTSS
ggplot(dat, aes(x=distanceToTSS, y=log2FoldChange)) +
  geom_point(alpha = 0.5) 


#########################################################################################
# Nearest Neighbor Analysis of DEG coding to non-coding genes
#########################################################################################

##-------------------------------------------------------## nearest distance of overlapping DEGs between Pol 2 and Pol 1

res.pol1.filtered$tx_id <- rownames(res.pol1.filtered)
res.pol2.filtered$tx_id <- rownames(res.pol2.filtered)
res.actd.filtered$tx_id <- rownames(res.actd.filtered)

pol1.df <- select(res.pol1.filtered, tx_id,seq_name, tx_seq_start, tx_seq_end, log2FoldChange, padj, tx_biotype, diffexpressed)
pol2.df <- select(res.pol2.filtered, tx_id,seq_name, tx_seq_start, tx_seq_end, log2FoldChange, padj, tx_biotype, diffexpressed)
actd.df <- select(res.actd.filtered, tx_id,seq_name, tx_seq_start, tx_seq_end, log2FoldChange, padj, tx_biotype, diffexpressed)

## filter counts on DEG data
plot <- merge(pol1.df, pol2.df, by = "tx_id")
plot <- merge(pol1.df, actd.df, by = "tx_id")
plot <- merge(pol2.df, actd.df, by = "tx_id")

## Set variables
lfc.up = 0.58
lfc.dn = -0.58
pvadj = 0.05
label.x = "Pol 1"
label.y = "Pol 2"

## label differentially expressed genes

plot$diffexpressed <- "NO"
plot$diffexpressed[plot$log2FoldChange.x > lfc.up & plot$log2FoldChange.y > lfc.up & plot$padj.x < pvadj & plot$padj.y < pvadj] <- "Upregulated"
plot$diffexpressed[plot$log2FoldChange.x < lfc.dn & plot$log2FoldChange.y < lfc.dn & plot$padj.x < pvadj & plot$padj.y < pvadj] <- "Downregulated"
plot$diffexpressed[plot$diffexpressed.x == "Downregulated" & plot$diffexpressed.y == "Upregulated"] <- paste0("Down in ",label.x, " | Up in ", label.y)
plot$diffexpressed[plot$diffexpressed.x == "Upregulated" & plot$diffexpressed.y == "Downregulated"] <- paste0("Up in ",label.x, " | Down in ", label.y)

#plot <- plot[plot$diffexpressed != "NO",] ## All Others down/upregulated

nndist.gene <- function(res, diffexp, biotype.1 = c("protein_coding"), biotype.2){
  
  require(plyranges)
  
  if ( diffexp == TRUE) {
    
    ## Get first results group
    res.1 <- res[res$diffexpressed == "Upregulated" | res$diffexpressed == "Downregulated" ,] 
    res.1 <- res.1[res.1$tx_biotype.x == biotype.1 ,]
    
    #res.1 <- res.1[res.1$diffexpressed == "Downregulated" ,] 
    #res.1 <- res.1[res.1$diffexpressed == "Upregulated" ,] 
    res.gr.1<- res.1 %>% 
      transform( seqnames= paste0("chr",res.1$seq_name.x), start = res.1$tx_seq_start.x, end = res.1$tx_seq_end.x)  %>% 
      as_granges()
    
    ## Get second results group 
    res.2 <- res[res$diffexpressed == "Upregulated" | res$diffexpressed == "Downregulated" ,] 
    #res.2 <- res.2[res.2$tx_biotype.x == biotype.2,]
    res.2 <- res.2[res.2$tx_biotype.x == "lncRNA" | res.2$tx_biotype.x == "retained_intron" | res.2$tx_biotype.x == "processed_pseudogene" | res.2$tx_biotype.x == "nonsense_mediated_decay"| res.2$tx_biotype.x == "protein_coding_CDS_not_defined",]
    
    #res.2 <- res.2[res.2$diffexpressed == "Downregulated" ,] 
    #res.2 <- res.2[res.2$diffexpressed == "Upregulated" ,] 
    res.gr.2<- res.2 %>% 
      transform( seqnames= paste0("chr",res.2$seq_name.x), start = res.2$tx_seq_start.x, end = res.2$tx_seq_end.x)  %>% 
      as_granges()
    
  } else { ## Get all genes
    
    res <- res[res$diffexpressed == "NO",] 
    
    ## Get first results group
    res.1 <- res[res$tx_biotype.x == biotype.1 ,]
    res.gr.1<- res.1 %>% 
      transform( seqnames= paste0("chr",res.1$seq_name.x), start = res.1$tx_seq_start.x, end = res.1$tx_seq_end.x)  %>% 
      as_granges()
    
    ## Get second results group 
    res.2 <- res[res$tx_biotype.x == biotype.2,]
    #res.2 <- res[res$tx_biotype.x == "lncRNA" | res$tx_biotype.x == "retained_intron" | res$tx_biotype.x == "processed_pseudogene" | res$tx_biotype.x == "nonsense_mediated_decay"| res$tx_biotype.x == "protein_coding_CDS_not_defined",]
    
    res.gr.2<- res.2 %>% 
      transform( seqnames= paste0("chr",res.2$seq_name.x), start = res.2$tx_seq_start.x, end = res.2$tx_seq_end.x)  %>% 
      as_granges()
  }
  
  x <- data.frame(plyranges::add_nearest_distance(res.gr.2 , res.gr.1, name = "distance"))
  
  return(x)
  
}

dge.nn.gene <- nndist.gene(plot, biotype.2="lncRNA", TRUE)
nn.gene <- nndist.gene(plot, biotype.2="lncRNA",FALSE)

nndist.gene <- function(res, diffexp, biotype.1 = c("lncRNA"), biotype.2){
  
  require(plyranges)
  
  if ( diffexp == TRUE) {
    
    ## Get Diff genes
    #res <- res[res$diffexpressed == "Upregulated" | res$diffexpressed == "Downregulated" ,] 
    
    ## Get first results group
    #res.1 <- res[res$gene_biotype == biotype.1 ,]
    res.1 <- res[res$tx_biotype == biotype.1 ,]
    
    #res <- res[res$diffexpressed == "Downregulated" ,] 
    #res <- res[res$diffexpressed == "Upregulated" ,] 
    res.gr.1<- res.1 %>% 
      #transform( seqnames= paste0("chr",res.1$seq_name), start = res.1$gene_seq_start, end = res.1$gene_seq_end)  %>% 
      transform( seqnames= paste0("chr",res.1$seq_name), start = res.1$tx_seq_start, end = res.1$tx_seq_end)  %>% 
      as_granges()
    
    ## Get second results group 
    res.2 <- res[res$gene_biotype == biotype.2,]
    #res.2 <- res[res$tx_biotype == biotype.2,]
    #res.2 <- res[res$tx_biotype == "lncRNA" | res$tx_biotype == "retained_intron" | res$tx_biotype == "processed_pseudogene" | res$tx_biotype == "nonsense_mediated_decay"| res$tx_biotype == "protein_coding_CDS_not_defined",]
   
    #res <- res[res$diffexpressed == "Downregulated" ,] 
    #res <- res[res$diffexpressed == "Upregulated" ,] 
    res.gr.2<- res.2 %>% 
      #transform( seqnames= paste0("chr",res.2$seq_name), start = res.2$gene_seq_start, end = res.2$gene_seq_end)  %>% 
      transform( seqnames= paste0("chr",res.2$seq_name), start = res.2$tx_seq_start, end = res.2$tx_seq_end)  %>% 
      as_granges()
    
  } else { ## Get all genes
    
    res <- res[res$diffexpressed != "Upregulated" & res$diffexpressed != "Downregulated" ,] 
    
    ## Get first results group
    #res.1 <- res[res$gene_biotype == biotype.1 ,]
    res.1 <- res[res$tx_biotype == biotype.1 ,]
    res.gr.1<- res.1 %>% 
      #transform( seqnames= paste0("chr",res.1$seq_name), start = res.1$gene_seq_start, end = res.1$gene_seq_end)  %>% 
      transform( seqnames= paste0("chr",res.1$seq_name), start = res.1$tx_seq_start, end = res.1$tx_seq_end)  %>% 
      as_granges()
    
    ## Get second results group 
    #res.2 <- res[res$gene_biotype == biotype.2,]
    res.2 <- res.pol2.filtered[res.pol2.filtered$tx_biotype == biotype.2,]
    #res.2 <- res[res$tx_biotype == biotype.2,]
    #res.2 <- res[res$tx_biotype == "lncRNA" | res$tx_biotype == "retained_intron" | res$tx_biotype == "processed_pseudogene" | res$tx_biotype == "nonsense_mediated_decay"| res$tx_biotype == "protein_coding_CDS_not_defined",]
    
    res.gr.2<- res.2 %>% 
      #transform( seqnames= paste0("chr",res.2$seq_name), start = res.2$gene_seq_start, end = res.2$gene_seq_end)  %>% 
      transform( seqnames= paste0("chr",res.2$seq_name), start = res.2$tx_seq_start, end = res.2$tx_seq_end)  %>% 
      as_granges()
  }
  
  x <- data.frame(plyranges::add_nearest_distance(res.gr.2 , res.gr.1, name = "distance"))
  
  return(x)
  
}

## Pol 1 
unique(res.pol1.filtered$tx_biotype)

dge.nn.gene <- nndist.gene(res.pol1.filtered, biotype.2="protein_coding", TRUE)
nn.gene <- nndist.gene(res.pol1.filtered, biotype.2="protein_coding",FALSE)

## Pol 2
dge.nn.gene <- nndist.gene(res.pol2.filtered,  biotype.2="lncRNA",TRUE)
nn.gene <- nndist.gene(res.pol2.filtered, biotype.2="lncRNA",FALSE)

## Act D
dge.nn.gene <- nndist.gene(res.actd.filtered,  TRUE)
nn.gene <- nndist.gene(res.actd.filtered, FALSE)

##-------------------------------------------------------## Plot dists

library(ggridges)

dge.nn.gene$cond <- rep("deg", nrow(dge.nn.gene))
nn.gene$cond <- rep("all", nrow(nn.gene))

data <- rbind(dge.nn.gene, nn.gene)
data$distance <- ifelse(data$distance ==0,log10(1), log10(data$distance))

dir = "/Users/lucascarter/Desktop/"
filename = "list.prot.cod"

write.csv(data, file=paste0(dir,filename, ".csv"), row.names=FALSE)

## Color scales
colors <- colorRampPalette(c("white","#798234","#D46780"))(8)
colors <- c("#D8DBC5", "#A5AB77", "#97794D", "#BA6E6A")

## Plot 
data %>% ## pol1.dist.deg.2.lnc.qtplot2  ## pol1.dist.lnc.2.deg.qtplot2
  ggplot(aes(x = distance, y = cond, fill = stat(quantile))) +
  stat_density_ridges(quantile_lines = F,calc_ecdf = TRUE,geom = "density_ridges_gradient", quantiles = 4) +theme_classic() +
  scale_y_discrete(breaks = c("deg","all" ), labels = c("DEGs", "non-DEGs"))+
  scale_x_continuous(breaks = seq(0, max(na.omit(data$distance))+ 1, 1))+
  scale_fill_manual(values=colors) + theme( legend.position="none",plot.title = element_text(size=12), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("Pol 1: distance from lncRNAs to Protein Coding genes") + xlab("Log10 Distance") + ylab("")


ggplot(data, aes(x = distance, y = cond)) + ## pol1.dist.deg.2.lnc.qtplot ## pol1.dist.lnc.2.deg.qtplot
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.75,
                      quantiles = 2)+theme_classic() + 
  scale_y_discrete(breaks = c("deg","all" ), labels = c("DEGs", "non-DEGs"))+
  scale_x_continuous(breaks = seq(0, max(na.omit(data$distance))+ 1, 1))+
  theme( legend.position="none",plot.title = element_text(size=12), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0))+
  ggtitle("Pol 1: distance from lncRNAs to Protein Coding genes") + xlab("Log10 Distance") + ylab("")


## get medians
mu <- ddply(data, "cond", summarise, mean=median(distance))

## ## pol1_2.dist.ncdeg.2.nonpc
ggplot(data, aes(x = distance, color = cond, fill = cond)) +
  geom_density(alpha = 0.8, color = NA) + 
  scale_color_manual(values = c("#A5AB77", "#BA6E6A"), aesthetics = c("color", "fill")) + theme_classic()+
  geom_vline(data=mu, aes(xintercept=mean, color=cond),linetype="dashed")+
  theme( legend.position="bottom",plot.title = element_text(size=10), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("Distance from DE protein coding genes to shared DE ncRNAs ") + xlab("Log10 Distance") + ylab("Density")


x = 119214625-200000
y = 119246827+200000

data <- res.pol1.filtered
dat <- data[data$seq_name == "8" & data$tx_seq_start >=  x  & data$tx_seq_end <= y,]

ggplot(dat, aes(x=tx_seq_start, y=log2FoldChange, color = tx_biotype, fill = tx_biotype)) +
  geom_point() 

#########################################################################################
# Nearest Neighbor Analysis of loops to non-coding genes
#########################################################################################

exps <- c("WT_HCT116_CTRL","Pol1_6hrs_Aux")

dir.loop <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/112123_HiC/112123_HiC/loop_analysis/Pol1_6hrs_Aux/mega/hiccups_results"
list.files(dir.loop)

##-------------------------------------------------------## Loop BED

dir.loop <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/112123_HiC/112123_HiC/loop_analysis/"
list.files(dir.loop)
exps[2]

## Convert results to Genomic Ranges data
loops <- as.data.frame(read.table(file.path(dir.loop,exps[2], '/mega/hiccups_results/',"merged_loops.bedpe"), sep= "\t", header= F))
#loops <- as.data.frame(read.table(file.path(dir.loop,exps[2], '/mega/hiccups_results/',"differential_loops1.bedpe"), sep= "\t", header= F)) ## Differential loop lists
loops <- loops %>% 
  transform( seqnames= paste0("chr",loops$V1), start = loops$V22, end = loops$V23)  %>% 
  #transform( seqnames= paste0("chr",loops$V1), start = loops$V2, end = loops$V6)  %>% 
  as_granges()

##-------------------------------------------------------## Loop anchors BED

dir.loop <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/112123_HiC/112123_HiC/loop_analysis/"
list.files(dir.loop)
exps[2]

## Convert results to Genomic Ranges data
loops.anch <- as.data.frame(read.table(file.path(dir.loop,exps[1], '/mega/hiccups_results/',"merged_loops.bedpe"), sep= "\t", header= F))

loops.5p <- loops.anch  %>% 
  transform( seqnames= paste0("chr",loops.anch$V1), start = loops.anch$V2, end = loops.anch$V3)  %>% 
  as_granges()

loops.3p <- loops.anch  %>% 
  transform( seqnames= paste0("chr",loops.anch$V1), start = loops.anch$V5, end = loops.anch$V6)  %>% 
  as_granges()

##-------------------------------------------------------## Get NN

nndist.loop <- function(res, loop.gr, diffexp, biotype){
  
  require(plyranges)
  
  if ( diffexp == TRUE) {
    
    ## Get Diff genes
    res <- res[res$diffexpressed == "Upregulated" | res$diffexpressed == "Downregulated" ,] 
    
    ## Get first results group
    #res <- res[res$gene_biotype == biotype ,]
    res <- res[res$tx_biotype == biotype ,]
    
    #res <- res[res$diffexpressed == "Downregulated" ,] 
    #res <- res[res$diffexpressed == "Upregulated" ,] 
    res.gr<- res %>% 
      #transform( seqnames= paste0("chr",res$seq_name), start = res$gene_seq_start, end = res$gene_seq_end)  %>% 
      transform( seqnames= paste0("chr",res$seq_name), start = res$tx_seq_start, end = res$tx_seq_end)  %>% 
      as_granges()
    
  } else { ## Get all genes
    
    res <- res[res$diffexpressed != "Upregulated" & res$diffexpressed != "Downregulated" ,] 
    
    ## Get first results group
    #res <- res[res$gene_biotype == biotype ,]
    res <- res[res$tx_biotype == biotype ,]
    res.gr<- res %>% 
      #transform( seqnames= paste0("chr",res$seq_name), start = res$gene_seq_start, end = res$gene_seq_end)  %>% 
      transform( seqnames= paste0("chr",res$seq_name), start = res$tx_seq_start, end = res$tx_seq_end)  %>% 
      as_granges()
    
  }
  
  x <- data.frame(plyranges::add_nearest_distance(res.gr , loop.gr, name = "distance"))
  
  return(x)
  
}

unique(res$gene_biotype)


## Pol 1 
dge.nn.loop <- nndist.loop(res.pol1.filtered, loop.gr = loops, TRUE, biotype= c("protein_coding"))
nn.loop <- nndist.loop(res.pol1.filtered, loop.gr = loops, FALSE, biotype= c("protein_coding"))

dge.nn.loop <- nndist.loop(res.pol1.filtered, loop.gr =loops, TRUE, biotype= c("lncRNA"))
nn.loop <- nndist.loop(res.pol1.filtered, loop.gr =loops, FALSE, biotype= c("lncRNA"))

## ActD 
dge.nn.loop <- nndist.loop(res.actd.filtered, loop.gr = loops.5p, TRUE, biotype= c("protein_coding"))
nn.loop <- nndist.loop(res.actd.filtered, loop.gr = loops.5p, FALSE, biotype= c("protein_coding"))

dge.nn.loop <- nndist.loop(res.actd.filtered, loop.gr =loops, TRUE, biotype= c("lncRNA"))
nn.loop <- nndist.loop(res.actd.filtered, loop.gr =loops, FALSE, biotype= c("lncRNA"))

## Pol2 
dge.nn.loop <- nndist.loop(res.pol2.filtered, loop.gr = loops, TRUE, biotype= c("protein_coding"))
nn.loop <- nndist.loop(res.pol2.filtered, loop.gr = loops, FALSE, biotype= c("protein_coding"))

dge.nn.loop <- nndist.loop(res.pol2.filtered, loop.gr =loops, TRUE, biotype= c("lncRNA"))
nn.loop <- nndist.loop(res.pol2.filtered, loop.gr =loops, FALSE, biotype= c("lncRNA"))

##-------------------------------------------------------## Plot dists

library(ggridges)

dge.nn.loop$cond <- rep("deg", nrow(dge.nn.loop))
nn.loop$cond <- rep("all", nrow(nn.loop))

data <- rbind(dge.nn.loop, nn.loop)
data$distance <- ifelse(data$distance ==0,log10(1), log10(data$distance))

## Color scales
colors <- colorRampPalette(c("white","#798234","#D46780"))(8)
colors <- c("#D8DBC5", "#A5AB77", "#97794D", "#BA6E6A")

## Plot 
data %>% ## pol1.dist.deg.2.lnc.qtplot2  ## pol1.dist.lnc.2.deg.qtplot2
  ggplot(aes(x = distance, y = cond, fill = stat(quantile))) +
  stat_density_ridges(quantile_lines = F,calc_ecdf = TRUE,geom = "density_ridges_gradient", quantiles = 4) +theme_classic() +
  scale_y_discrete(breaks = c("deg","all" ), labels = c("DEGs", "non-DEGs"))+
  scale_x_continuous(breaks = seq(0, max(na.omit(data$distance))+ 1, 1))+
  scale_fill_manual(values=colors) + theme( legend.position="none",plot.title = element_text(size=12), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("distance from loops to protein coding genes") + xlab("Log10 Distance") + ylab("")


ggplot(data, aes(x = distance, y = cond)) + ## pol1.dist.deg.2.lnc.qtplot ## pol1.dist.lnc.2.deg.qtplot
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.75,
                      quantiles = 2)+theme_classic() + 
  scale_y_discrete(breaks = c("deg","all" ), labels = c("DEGs", "non-DEGs"))+
  scale_x_continuous(breaks = seq(0, max(na.omit(data$distance))+ 1, 1))+
  theme( legend.position="none",plot.title = element_text(size=12), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0))+
  ggtitle("Pol 1: distance from lncRNAs to Loops") + xlab("Log10 Distance") + ylab("")


ggplot(data, aes(x=distance, y=log2FoldChange, color = cond, fill = cond)) +
  geom_point() +
  scale_color_manual(values = c("#A5AB77", "#BA6E6A"), aesthetics = c("color", "fill"))

#####################################################################################
# Intronic coverage plotting                         
#####################################################################################
library(dplyr)
library(tidyr)
library(pheatmap)
library(gridExtra)
library(ggplotify)

cvg <- readRDS('/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Intron-reads/intron_coverage/sumcvg.rds')

pivot_data <- function(size, feature, name, conds) {

    df <- data.frame()
    for (i in 1:length(conds)){
    
    cond <- conds[i]
    
    group <- cvg %>% dplyr::filter(Group == cond)
    
    tbl<- group %>% dplyr::filter(cat == size & feature_type == feature)  %>% select(-Group, -cat, -feature_type, -bumpy) %>% pivot_wider(names_from = bin, values_from = lumpy)
    tbl <- data.frame(tbl)
    rownames(tbl) <- cond
    
    
    df <- rbind(df, tbl)
    }

colnames(df) <- seq(1,20)

return(df)

}

shr.int <- pivot_data(size = "short", feature = "intron", conds = dput(unique(cvg$Group)))
med.int <- pivot_data(size = "regular", feature = "intron", conds = dput(unique(cvg$Group)))
lng.int <- pivot_data(size = "long", feature = "intron", conds = dput(unique(cvg$Group)))
shr.ex <- pivot_data(size = "short", feature = "exon", conds = dput(unique(cvg$Group)))
med.ex <- pivot_data(size = "regular", feature = "exon", conds = dput(unique(cvg$Group)))
lng.ex <- pivot_data(size = "long", feature = "exon", conds = dput(unique(cvg$Group)))

colors <- colorRampPalette(c("#798234","white","#D46780"))(64)

mats <- list(shr.int, med.int, lng.int, shr.ex, med.ex, lng.ex)
p <- list()
for(i in 1:length(mats)){
  
  mat <- mats[[i]]
  print(mat)
  p[[i]] <- as.ggplot(pheatmap(mat , display_numbers = F, number_color = "black", color = colors,
                     cellwidth = 8, cellheight = 8,labels_row = c("ActD 6hrs", "CTRL DMSO", "Pol1 8hrs", "Pol2 8hrs"),
                     fontsize = 8,  cluster_rows = F,  cluster_cols = F))
 
}

grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],
             nrow=2, ncol = 3,  top = "Intronic/Exonic Coverage",vp=viewport(width=1, height=0.5, clip = TRUE)
)

library(cowplot)
plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]])

## Plot panel of all variables
cvg <- na.omit(cvg)
colors <- c("#D8DBC5", "#A5AB77", "#97794D", "#BA6E6A")
colors <- c("#009392","#cf597e","#e88471","#39b185")
ggplot(cvg, aes(x = bin, y = lumpy, colour = Group, group = Group)) +
  geom_line(size = 1) +
  facet_grid(feature_type ~ cat, scales = "free_y") + ## panel of conditions
  scale_x_continuous(limits = c(1,20)) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size=16), 
        text = element_text(size=16, family="Arial")
  ) +
  labs(y = "Normalized Log-Coverage", 
       x = "Binned Gene Body (5' -> 3')")


## Plot panel of aggregate values
cvg2 <- cvg %>% select(-cat) %>% group_by(feature_type,bin, Group) %>% summarise( sum = mean(lumpy))

ggplot(cvg2, aes(x = bin, y = sum, colour = Group, group = Group)) +
  geom_line(size = 1) +
  facet_grid(rows = vars(feature_type),scales = "free_y") + ## panel of conditions
  scale_x_continuous(limits = c(1,20)) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size=16), 
        text = element_text(size=16, family="Arial")
  ) +
  labs(y = "Normalized Log-Coverage", 
       x = "Binned Gene Body (5' -> 3')")


#####################################################################################
#  check loops on introns                      
#####################################################################################

getOverlaps <- function(exps, dir.loop, dir.int, interval, cond, anchors, degs){
  
  ## exps: vector of experiment labels # exps <- c("WT_HCT116_CTRL","Pol1_6hrs_Aux","Pol2_6hrs_Aux")
  ## dir.loop: root directory of loops # dir.loop <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/112123_HiC/112123_HiC/loop_analysis/"
  ## dir.int: root director of genomic ranges intervals # dir.int <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations/Intron_GTFs"
  ## cond: a string for the condition name
  ## anchors: boolean for if using 5' and 3' loop anchors instead of full loop
  ## degs: boolean if using DEG list. values set to 0.05 and 0.58
  ## Loop file must be in HICCUPS format
  
  require(plyranges)
  require(dplyr)

##-------------------------------------------------------## Loop BED file
loops <- as.data.frame(read.table(file.path(dir.loop,exps[1], '/mega/hiccups_results/',"merged_loops.bedpe"), sep= "\t", header= F))

##-------------------------------------------------------## Loop anchors overlap
if ( anchors == T ){
  
loops.5p <- loops %>% transform( seqnames= paste0("chr",V1), start = V2, end = V3)  %>% as_granges()
loops.3p <- loops  %>% transform( seqnames= paste0("chr",V1), start = V5, end = V6)  %>% as_granges()

ints <- as.data.frame(read.table(file.path(dir.int,interval), sep= "\t", header= T))

int.overlaps.3p <- data.frame(join_overlap_inner(loops.3p, ints))
int.overlaps.5p <- data.frame(join_overlap_inner(loops.5p, ints))

int.out <- rbind(int.overlaps.3p, int.overlaps.5p)
int.out <- int.out  %>% dplyr::mutate(dups = paste0(V1, "_", V22, "_",V23)) %>% dplyr::distinct(dups, .keep_all = T) ## Filter dups
int.out$cond <- rep(cond, nrow(int.out))

return(int.out)
##-------------------------------------------------------## loop DEG overlaps
} else if ( degs == T ) {
  
loops <- loops %>% transform( seqnames= paste0("chr",V1), start = V22, end = V23) %>% as_granges()

ints <- as.data.frame(read.table(file.path(dir.int,interval), sep= "\t", header= T))
ints <- ints %>% filter(abs(logFC) > 0.58, ints$adj.P.Val < 0.05) ## filter ints for DE
#ints <- ints %>% filter(logFC < 0.58, ints$adj.P.Val < 0.05) ## filter ints for DE
#ints <- ints %>% filter(logFC > 0.58, ints$adj.P.Val < 0.05) ## filter ints for DE

ints <- ints  %>% transform( seqnames= paste0(Chr), start = Start , end = End)  %>% as_granges()

int.out <- data.frame(join_overlap_inner(loops, ints))
int.out <- int.out  %>% dplyr::mutate(dups = paste0(V1, "_", V22, "_",V23)) %>% dplyr::distinct(dups, .keep_all = T) ## Filter dups
int.out$cond <- rep(cond, nrow(int.out))

return(int.out)
##-------------------------------------------------------## loop interval overlaps
} else {

loops <- loops %>% transform( seqnames= paste0("chr",V1), start = V22, end = V23)  %>% as_granges()

ints <- as.data.frame(read.table(file.path(dir.int,interval), sep= "\t", header= T))
ints <- ints  %>% transform( seqnames= paste0(Chr), start = Start , end = End)  %>% as_granges()

int.out <- data.frame(join_overlap_inner(loops, ints))
int.out <- int.out  %>% dplyr::mutate(dups = paste0(V1, "_", V22, "_",V23)) %>% dplyr::distinct(dups, .keep_all = T) ## Filter dups
int.out$cond <- rep(cond, nrow(int.out))

return(int.out)
  
}

}

## define repeat variables
dir.loop <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/112123_HiC/112123_HiC/loop_analysis/"
dir.int <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations/Intron_GTFs"
dir.deg <-"/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Intron-reads/Analyses/Results"
exps <- c("WT_HCT116_CTRL","Pol1_6hrs_Aux","Pol2_6hrs_Aux")

## call function
ints.wt.all <- getOverlaps(exps = exps[1], dir.loop = dir.loop, dir.int=dir.int, interval="intron.txt", cond = "WT",anchors = F, degs = F)
ints.p1.all <- getOverlaps(exps = exps[2], dir.loop = dir.loop, dir.int=dir.int, interval="intron.txt", cond = "Pol 1", anchors = F, degs = F)

ints.wt <- getOverlaps(exps = exps[1], dir.loop = dir.loop, dir.int=dir.deg, interval="pol1.intron.dge.txt", cond = "WT", anchors = F, degs = T)
ints.p1 <- getOverlaps(exps = exps[2], dir.loop = dir.loop, dir.int=dir.deg, interval="pol1.intron.dge.txt", cond = "Pol 1", anchors = F, degs = T)

ints.wt <- getOverlaps(exps = exps[1], dir.loop = dir.loop, dir.int=dir.deg, interval="pol2.intron.dge.txt", cond = "WT",anchors = F, degs = T)
ints.p2 <- getOverlaps(exps = exps[3], dir.loop = dir.loop, dir.int=dir.deg, interval="pol2.intron.dge.txt", cond = "Pol 2", anchors = F, degs = T)

ints.p1.dn <- getOverlaps(exps = exps[2], dir.loop = dir.loop, dir.int=dir.deg, interval="pol1.intron.dge.txt", cond = "Pol 1 Down", anchors = F, degs = T)
ints.p1.up <- getOverlaps(exps = exps[2], dir.loop = dir.loop, dir.int=dir.deg, interval="pol1.intron.dge.txt", cond = "Pol 1 UP", anchors = F, degs = T)


##-------------------------------------------------------## Plot

int.plot <- rbind(ints.wt, ints.p1)
int.plot <- rbind(ints.wt.all , ints.p1.all)
int.plot <- rbind(ints.p1.dn , ints.p1.up)
int.plot <- rbind(ints.wt, ints.p2)
int.plot.f <- data.frame(chrom = int.plot$seqnames,count = rep(1, nrow(int.plot)),cond = int.plot$cond, dist = abs(int.plot$V2 - int.plot$V6), obs = int.plot$V12, lfc = int.plot$logFC)
max(int.plot$logFC)
library(ggplot2)
library(ggpubr)
library(forcats)

c("#5F4690","#1D6996","#38A6A5","#0F8554","#73AF48","#EDAD08","#E17C05","#CC503E","#94346E","#6F4070","#994E95","#666666","#009392","#cf597e","#e88471","#39b185")

pos <- position_jitter(height = 0.05, width = 0.05, seed = 1) 
## Plot by quartile # pol1.lfc.vs.pol2.lfc.allgenes
ggplot(int.plot.f, aes(x=log10(abs(lfc)), y=log10(dist), color = cond)) +
  geom_point(position = pos, alpha=0.5, size=1)+ theme_pubr()+
  scale_color_manual(values=c("#EDAD08","#5F4690"), aesthetics = c("color", "fill"))+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  xlab("Log10 abs LFC")+ ylab("Log10 Observed")

int.plot$V1 <- ifelse(int.plot$V1 == "X", as.numeric(23), as.numeric(int.plot$V1))
cats <- data.frame(chrom = int.plot$seqnames,factor = as.numeric(int.plot$V1), count = rep(1, nrow(int.plot)),cond = int.plot$cond, dist = abs(int.plot$V2 - int.plot$V6), obs = int.plot$V12)
tbl <- cats  %>% group_by(chrom,cond) %>% summarise( sum = sum(count)) %>% ungroup()
tbl$factor<- as.numeric(ifelse(tbl$chrom == "chrX" , 23, sapply(strsplit(as.character(tbl$chrom), "chr"), `[`, 2)))

tbl %>%
  mutate(chrom= fct_reorder(as.character(chrom),as.numeric(factor))) %>% 
  ggplot(aes(y=sum, x=chrom, fill=cond)) + theme_classic() +
  geom_bar(stat="identity", position="dodge",color= "black") + scale_fill_manual(values = c( "white", "black"))+
  xlab("Chroms")+ ylab("Number of DE Loops")

#########################################################################################
# POLR3A motif results
#########################################################################################
library(plyranges)
library(forcats)

dir.int <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations/Intron_GTFs"
dir.deg <-"/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Intron-reads/Analyses/Results"
dir.mot <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/pol3_consensus"
dir.mot2 <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/pol3_consensus/meme_chip_out/fimo_out_1"
list.files(dir.mot2)

## get int ranges
ints.p1 <- as.data.frame(read.table(file.path(dir.deg,"pol1.intron.dge.txt"), sep= "\t", header= T)) %>% filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% mutate(count = rep(1, 1219)) %>% transform( seqnames= paste0(Chr), start = Start, end = End)  %>% as_granges()
ints.p2 <- as.data.frame(read.table(file.path(dir.deg,"pol2.intron.dge.txt"), sep= "\t", header= T)) %>% filter(abs(logFC) > 0.58, adj.P.Val < 0.05)  %>% mutate(count = rep(1, 2627)) %>% transform( seqnames= paste0(Chr), start = Start, end = End)  %>% as_granges()
ints <- as.data.frame(read.table(file.path(dir.int,"intron.txt"), sep= "\t", header= T))  %>% mutate(count = rep(1, 296770)) %>% transform( seqnames= paste0(Chr), start = Start, end = End)  %>% as_granges()

## Get homer motif file
mots <- as.data.frame(read.table(file.path(dir.mot,"POLR3A.sites.5pmasked.intersect.hg38.bed"), sep= "\t", header= F)) %>% transform( seqnames= paste0(V1), start = V2, end = V3)  %>% as_granges()
mots <- as.data.frame(read.table(file.path(dir.mot2,"fimo.gff"), sep= "\t", header= F)) %>% transform( seqnames= paste0(V1), start = V4, end = V5)  %>% as_granges()

## Get motif overlaps for each group
ints.p1.ov <- data.frame(join_overlap_inner(ints.p1, mots)) 
ints.p1.ov <- ints.p1.ov %>% mutate(count = rep(1, nrow(ints.p1.ov)))

ints.p2.ov <- data.frame(join_overlap_inner(ints.p2, mots))
ints.p2.ov <- ints.p2.ov %>% mutate(count = rep(1, nrow(ints.p2.ov)))

ints.all.ov<- data.frame(join_overlap_inner(ints, mots))
ints.all.ov<- ints.all.ov%>% mutate(count = rep(1, nrow(ints.all.ov)))

## Build table control tables
p1.ints.tbl <- data.frame(ints.p1) %>% group_by(Chr) %>% summarise( sum = sum(count))
p2.ints.tbl <- data.frame(ints.p2) %>% group_by(Chr) %>% summarise( sum = sum(count))
ints.tbl <- data.frame(ints) %>% group_by(Chr) %>% summarise( sum = sum(count))

## Build treatment tables
p2.tbl <- ints.p2.ov  %>% group_by(Chr) %>% summarise( sum = sum(count)) 
p1.tbl <- ints.p1.ov  %>% group_by(Chr) %>% summarise( sum = sum(count)) 
all.tbl <- ints.all.ov  %>% group_by(Chr) %>% summarise( sum = sum(count)) 

## ratios
p1 <- data.frame(ratio = p1.tbl$sum/p1.ints.tbl$sum[1:23], chrom = p1.tbl$Chr, cond  = rep("Pol1 DE Introns", nrow(p1.tbl)) )
p2 <- data.frame(ratio = p2.tbl$sum/p2.ints.tbl$sum[1:23], chrom = p2.tbl$Chr, cond  = rep("Pol2 DE Introns", nrow(p2.tbl)) )
all <- data.frame(ratio = all.tbl$sum/ints.tbl$sum, chrom = all.tbl$Chr, cond  = rep("All Introns", nrow(all.tbl)) )
all <- all[c(1:23),]

plot <- rbind(p1,p2,all)
plot <- mutate(plot, factor= as.numeric(ifelse(chrom == "chrX" , 23, sapply(strsplit(as.character(chrom), "chr"), `[`, 2))))

c("#5F4690","#1D6996","#38A6A5","#0F8554","#73AF48","#EDAD08","#E17C05","#CC503E","#94346E","#6F4070","#994E95","#666666","#009392","#cf597e","#e88471","#39b185")

## Plot loops for each chromosome
plot %>% 
  mutate(chrom = fct_reorder(chrom, as.numeric(factor))) %>% 
  ggplot( aes(fill=cond, y=ratio, x=reorder(chrom, factor), group = cond)) + 
  geom_bar(position=position_dodge(), stat="identity") + theme_classic()+
  guides(x=guide_axis(angle = 45)) + 
  scale_y_continuous(labels = scales::percent) + coord_cartesian(ylim = c(0, 1))+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(size = 8)) +
  scale_fill_manual(values = c("#5F4690","#38A6A5","#EDAD08"))+ ggtitle("POLR3A Binding Profile")+
  labs(x = 'Chromosomes', y = '% of Total')



#########################################################################################
# Generate heatmap of top Differential results 
#########################################################################################

#https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/
## As gene.n (number of tiles) increases, tile.width must increase, to keep image proportional
Heat_Map_Counts <- function(results, title, subtitle, caption, sample.n = 1:8, padj = 0.05, lfc.abs = 1, gene.n=gene.n, 
                            n_colors=n_colors, breaks=breaks, labels = as.character(labels), tile.width=tile.width,
                            axis.text.x = element_text( color="black", size=7, family = "Arial Narrow", angle = 90, vjust = 0.5, hjust=1),
                            axis.text.y = axis.text.y) {
  
  require(ggplot2) # ggplot() for plotting
  require(dplyr) # data reformatting
  require(tidyr) # data reformatting
  require(stringr) # string manipulation
  require(viridis) # colors package
  require(DESeq2) # DESeq2
  
  
  ## Regularized log transformation, count extraction, and normalization
  rld <- DESeq2::rlog(dds, blind = FALSE)
  assay <-t(scale(t(assay(rld)[,sample.n])))
  
  ## Remove .## version of gene in ensembl ID
  rownames(assay) <- gsub("\\..*","",rownames(assay))
  
  ## Organize genes by significance and extract
  sig.genes <- rownames(results[results$padj <= padj & abs(results$log2FoldChange) > abs(lfc.abs),])
  results <- as.data.frame(results[rownames(results) %in% sig.genes,])
  results <- results %>% 
    dplyr::distinct(symbol, .keep_all = T)
  
  ## Order genes by LFC, pAdj, and extract N genes for heatmap
  #results <- results[order(results$padj, decreasing = F),]
  #results <- as.data.frame(results[order(results$log2FoldChange, decreasing = FALSE),])[gene.n,] # order the results DOWNREG
  #results <- as.data.frame(results[order(results$log2FoldChange, decreasing = TRUE),])[gene.n,] # order the results UPREG
  results <- as.data.frame(results)[gene.n,] # unordered results for clustering
  
  ## Prepare data frame for heatmapping 
  top.heat <- as.data.frame(assay[rownames(assay)%in%rownames(results),]) %>% 
    dplyr::mutate(label = rownames(results)) %>%
    dplyr::mutate(symbol = results$symbol) %>% # make rownames a new column
    tidyr::gather(key="samples",value="value", -symbol, -label) %>% # convert data to long format
    stats::setNames(c("label", "symbol","samples","value")) %>% # rename columns
    dplyr::mutate(regions=factor(symbol)) %>% # convert factors
    dplyr::mutate(samples=factor(samples)) %>%  # convert factors
    dplyr::mutate(value=as.numeric((value))) #convert value to numeric (also converts '-' to NA, gives a warning)
  
  ## Compute a distance calculation on both dimensions of the matrix for clustering
  #distance.gene <- dist(as.data.frame(assay[rownames(results),])) ## average distance method
  distance.gene <- as.dist((1 - cor(t(as.data.frame(assay[rownames(results),]))))/2) ## correlation method
  distance.sample <- dist(t(as.data.frame(assay[rownames(results),])))
  
  # Cluster based on the distance calculations
  #cluster.gene <- hclust(distance.gene, method="average") ## average distance method
  cluster.gene=hclust(distance.gene) ## correlation method
  cluster.sample <- hclust(distance.sample, method="average")
  
  # Re-factor sample and genes for ggplot2
  top.heat$label <- factor(top.heat$label, levels=cluster.gene$labels[cluster.gene$order])
  top.heat$samples <- factor(top.heat$samples, levels=cluster.sample$labels[cluster.sample$order])
  top.heat$regions <- factor(top.heat$regions, levels=top.heat$regions[cluster.gene$order])
  
  # upregulated color: #D22B2B
  # downregulated color: #0047AB
  ## Apply viridis_pal function to generate a palette
  #palette <- colorRampPalette(c("#D22B2B", "#0047AB"))(n_colors)
  palette <- colorRampPalette(c("#0047AB", "white", "#D22B2B"))(n_colors)
  
  ## Bin heat map values using cut() and factor
  top.heat <- top.heat %>% # create a new variable from count
    mutate(countfactor=cut(value,breaks=c(breaks), labels=c(labels))) %>% # change level order
    mutate(countfactor=factor(as.character(countfactor),levels=rev(levels(countfactor)))) %>% # factor bins
    arrange(factor(top.heat$label, levels=cluster.gene$labels[cluster.gene$order])) %>% # factor labels
    arrange(factor(top.heat$regions, levels=top.heat$regions[cluster.gene$order])) # factor regions
  
  
  ## Plot heat map of top 1000
  heat.map <- top.heat %>% 
    ggplot(aes(x=samples,y=regions,fill=countfactor))+
    geom_tile(size=1)+coord_fixed(ratio=tile.width)+
    scale_fill_manual(values=c(palette),na.value = "grey90")+ ## color scale 2
    labs(title=title, subtitle=subtitle, caption=caption) +
    scale_y_discrete(name=NULL)+
    scale_x_discrete(name=NULL)+
    theme(axis.text.y = axis.text.y ,
          axis.text.x = axis.text.x ,
          plot.title = element_text(hjust = 0.5, color="black", face = "bold", size=15,
                                    family = "Arial Rounded MT Bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  
  return(heat.map)
  
}

## Call function for N=1000 and N=100
Heat_Map_Counts(res.pol1, "Top 100 Differentially Expressed Genes", "", "", gene.n=1:100, tile.width = 0.4, breaks=seq(from=-1.75, to=1.75, by=0.25),
                n_colors=14, labels=seq(from=-1.75, to=1.5, by=0.25), axis.text.y = element_text( color="black", size=4.5, family = "Arial Narrow"))
Heat_Map_Counts(res.pol1, "Top 1000 Differentially Expressed Genes", "", "", gene.n=1:673, tile.width = .04, breaks=seq(from=-2, to=1.75, by=0.25),
                n_colors=15, labels=seq(from=-2, to=1.5, by=0.25),axis.text.y =element_blank())
Heat_Map_Counts(res.12hr, "All Differentially Expressed Genes", "", "", gene.n=1:2039, tile.width = .01, breaks=seq(from=-2.25, to=2.5, by=0.25),
                n_colors=18, labels=seq(from=-2.25, to=2.25, by=0.25),axis.text.y =element_blank())

#####################################################################################
# Plot tracks                                     
#####################################################################################

library(biomaRt)
library(Gviz)
library(plyranges)
library(dplyr)
library(rtracklayer)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

## BED files
ctcf <- read.delim('/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations/ctcf/ENCFF463FGL.bed',h = F, sep = '\t', stringsAsFactors = F)
rad21 <- read.delim('/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations/rad21/ENCFF391AAM.bed',h = F, sep = '\t', stringsAsFactors = F)
pol2ps5 <- read.delim('/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations/polr2a-ps5/ENCFF229YMU.bed',h = F, sep = '\t', stringsAsFactors = F)

cov <- import.bw("/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Backman19_10.18.2023/totalRNA_BAMs/pol2.merge.bigWig",as="GRanges") 
## BIGWIGS
k4me3.bw <- '/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations/H3K4me3/ENCFF394IMO.bigWig'
k9me3.bw <- '/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations/H3K9me3/ENCFF572IBD.bigWig'

## RNAseq
actd.rna <- '/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/bigwigs/actd.merge.rpgc.bigWig'
wt.rna <- '/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/bigwigs/wt.merge.bigWig'
pol1.rna <- import.bw('/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/bigwigs/pol1.merge.rpgc.bigWig',as="GRanges")
pol2.rna <- import.bw('/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/bigwigs/pol2.merge.rpgc.bigWig',as="GRanges") 

## Retrieve genes

values=c("ENSG00000147676") ## Mal2

gb <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','chromosome_name', "start_position","end_position","exon_chrom_start","exon_chrom_end", "strand"), filters =
              "ensembl_gene_id", values=values, mart=ensembl)

colnames(gb) <- c('ENSG','ENST','chrom', 'txStart' , 'txEnd','exonStart' , 'exonEnd', 'strand')

gb$chrom = gsub(gb$chrom, pattern = '^', replacement = 'chr')
gb$strand = ifelse(gb$strand == 1, '+',"-")

## Plot gene track 
biomTrack <- BiomartGeneRegionTrack(genome = "hg38",
                                    gene=values, ## whole myosin chain
                                    name = "ENSEMBL", biomart = ensembl, fill = "black", color="black", cex.feature = 0.75,
                                    fontsize=12,fontcolor.title="white",background.title = "black", fontcolor.group="black")

##-- Plotting bigwig tracks --##
chrom = 'chr8'
## Colors
# green "#798234"
# K9 red "#cf597e"
# K4 blue "#3969AC"
# K27 orange "#E68310"
##

dTrack1 <- DataTrack(range = k9me3.bw, genome = "hg38",
                     chromosome = chrom, 
                     name = "H3K9me3", type = "hist", 
                     #window = -1, windowSize = 10, 
                     ylim=c(2,10),col.histogram="#cf597e", background.title = "#cf597e")

dTrack2 <- DataTrack(range = k4me3.bw, genome = "hg38",  
                     chromosome = chrom, 
                     name = "H3K4me3", type = "hist",
                     #window = -1, windowSize = 10, 
                     ylim=c(2,10), col.histogram="#E68310", background.title = "#E68310")

dTrack3 <- DataTrack(range = pol2ps5, genome = "hg38",
                     chromosome = chrom,
                     name = "Pol2-PS5",
                     #window = -1, windowSize = 10, 
                     ylim=c(2,30),col.histogram="#3969AC", background.title = "#3969AC")

dTrack4 <- DataTrack(range = rad21, genome = "hg38",
                     chromosome = chrom, 
                     name = "Rad21", type = "line", 
                     #window = -1, windowSize = 500,
                     ylim=c(0,5),col.histogram="black", background.title = "black")

dTrack5 <- DataTrack(range = ctcf, genome = "hg38",
                     chromosome = chrom, 
                     name = "CTCF", type = "line", 
                     #window = -1, windowSize = 5000, ylim=c(0,8),
                     col.histogram="black", background.title = "black")

dTrack6 <- DataTrack(range = pol1.rna, genome = "hg38",
                     chromosome = chrom, 
                     name = "Pol 1 Transcription", type = "hist", 
                     window = -1, windowSize = 500,
                     ylim=c(0,5),col.histogram="black", background.title = "black")

dTrack7 <- DataTrack(range = pol2.rna, genome = "hg38",
                     chromosome = chrom, 
                     name = "Pol 2 Transcription", type = "hist", 
                     window = -1, windowSize = 5000, ylim=c(0,8),
                     col.histogram="black", background.title = "black")

dTrack8 <- DataTrack(range = actd.rna, genome = "hg38",
                     chromosome = chrom, 
                     name = "ActD Transcription", type = "hist", 
                     window = -1, windowSize = 500,
                     ylim=c(0,5),col.histogram="black", background.title = "black")

dTrack9 <- DataTrack(range = wt.rna, genome = "hg38",
                     chromosome = chrom, type = "histogram", importFunction=Gviz:::.import.bw,#window = -1, windowSize = 10,
                     col.histogram="black", background.title = "black")

x = 119165034	-100000

y = 119245673 + 100000

## Plot full locus
plotTracks(list(biomTrack, dTrack1, dTrack2, dTrack3, dTrack4, dTrack5, dTrack6, dTrack7, dTrack8, dTrack9), from = x, to =y,  collapseTranscripts="meta", showID=T, transcriptAnnotation="symbol")
plotTracks(list(biomTrack, dTrack1, dTrack2, dTrack6, dTrack9), from = x, to =y,  collapseTranscripts="meta", showID=T, transcriptAnnotation="symbol")
plotTracks(list( dTrack9), from = x, to =y, genome="hg38")
#####################################################################################
# DE gene pheatmap                                         
#####################################################################################

# upregulated color: #D22B2B
# downregulated color: #0047AB

## Set colors for annotations
q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = "E") # E= Civis
v_colors

Gene_pHeat_Map <- function(results, title, gene.n = gene.n, padj = 0.05, lfc.abs = 0.58, sample.n = 1:8, c.width = c.width, c.height = c.height, f.size.row=4,
                           f.size.col=10, f.size=10, cols=colorRampPalette(c("#0047AB", "white", "#D22B2B"))(30), 
                           meta.annotation=c("Conditions", "Replicate"),
                           anno.colors = list(Replicate = c("1" = "#000004FF", "2" = "#5E626EFF"), 
                                             Conditions = c("0hr" = "#00204DFF", "12hr" = "#00336FFF", "48hr" = "#39486BFF", "144hr" = "#575C6DFF"))
                           ) {
  require(pheatmap)
  require(viridis)
  require(dendsort)
  require(DESeq2)
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  ## Remove .## version of gene in ensembl ID
  #rownames(results) <- gsub("\\..*","",rownames(results))
  #rownames(dds) <- gsub("\\..*","",rownames(dds))
  
  
  ## Regularized log transformation, count extraction, and normalization
  rld <- DESeq2::rlog(dds, blind = FALSE)
  assay <-t(scale(t(assay(rld)[,sample.n])))
  
  ## Remove .## version of gene in ensembl ID
  rownames(assay) <- gsub("\\..*","",rownames(assay))
  
  ## Organize genes by significance and extract
  sig.genes <- rownames(results[results$padj <= padj & abs(results$log2FoldChange) > abs(lfc.abs),])
  results <- as.data.frame(results[rownames(results) %in% sig.genes,])
  results <- results %>% 
    dplyr::distinct(symbol, .keep_all = T)
  
  ## Order genes by  pAdj, and extract N genes for heatmap
  results <- results[order(results$padj, decreasing = F),]
  print(results)
  #results <- as.data.frame(results[order(results$log2FoldChange, decreasing = FALSE),])[gene.n,]
  results <- as.data.frame(results)[gene.n,]
  results
  
  ## Prepare data frame for heatmapping 
  top.heat <- as.data.frame(assay[rownames(assay)%in%rownames(results),]) 
  top.heat
  ## Compute a distance calculation on both dimensions of the matrix for clustering
  #distance.gene <- dist(as.data.frame(assay[rownames(results),])) ## average distance method
  distance.gene <- as.dist((1 - cor(t(as.data.frame(assay[rownames(results),]))))/2) ## correlation method
  distance.sample <- dist(t(as.data.frame(assay[rownames(results),])))
  
  # Cluster based on the distance calculations
  #cluster.gene <- hclust(distance.gene, method="average") ## average distance method
  cluster.gene=hclust(distance.gene) ## correlation method
  cluster.sample <- hclust(distance.sample, method="average")
  
  ## Annotation information in data frame object
  dfcol <- as.data.frame(colData(dds)[,meta.annotation])

  # plot heat map
  p<-pheatmap(top.heat, 
              color = cols,
              fontsize = f.size,
              fontsize_row = f.size.row, 
              fontsize_col = f.size.col,
              cellwidth = c.width,
              cellheight = c.height,
              cluster_cols= F,
              #cluster_cols= cluster.sample,
              cluster_rows= cluster.gene,
              main = title,
              annotation_colors = anno.colors,
              labels_row = results$symbol,
              show_rownames = T,
              border_color = NA,
              annotation_legend = TRUE,
              annotation_col = dfcol
              )
 
  return(p)
  
}

Gene_pHeat_Map(res.12hr.ov, "", gene.n=1:100, c.width = 7, c.height = 4) #12hr_inLADS_heatmap_top100
Gene_pHeat_Map(res.12hr.non, "", gene.n=1:100, c.width = 7, c.height = 4) #12hr_outsideLADS_heatmap_top100

## Use these for new plots
Gene_pHeat_Map(res.48hr.ov, "LAD overlapping", gene.n=1:100, c.width = 9, c.height = 6) #48hr_inLADS_heatmap_top100
Gene_pHeat_Map(res.48hr.non, "LAD non-overlapping", gene.n=1:100, c.width = 9, c.height = 6) #48hr_outsideLADS_heatmap_top100
##

Gene_pHeat_Map(res.144hr.ov, "", gene.n=1:100, c.width = 7, c.height = 4) #144hr_inLADS_heatmap_top100
Gene_pHeat_Map(res.144hr.non, "", gene.n=1:100, c.width = 7, c.height = 4) #144hr_outsideLADS_heatmap_top100

res <- res.144hr.non
res <-res%>%
  filter(abs(log2FoldChange) > 0.58 )%>%
  filter(padj < 0.05 )

Gene_pHeat_Map(res.12hr.ov, "", gene.n=1:500, c.width = 7, c.height = 0.75) #12hr_inLADS_heatmap_top500
Gene_pHeat_Map(res.12hr.non, "", gene.n=1:500, c.width = 7, c.height = 0.75) #12hr_outsideLADS_heatmap_top500

Gene_pHeat_Map(res.48hr.ov, "", gene.n=1:500, c.width = 7, c.height = 0.75) #48hr_inLADS_heatmap_top500
Gene_pHeat_Map(res.48hr.non, "", gene.n=1:500, c.width = 7, c.height = 0.75) #48hr_outsideLADS_heatmap_top500

Gene_pHeat_Map(res.144hr.ov, "", gene.n=1:400, c.width = 7, c.height = 0.75) #144hr_inLADS_heatmap_top400
Gene_pHeat_Map(res.144hr.non, "", gene.n=1:500, c.width = 7, c.height = 0.75) #144hr_outsideLADS_heatmap_top500

rm(assay)

#####################################################################################
# DE gene pheatmap                                         
#####################################################################################

# upregulated color: #D22B2B
# downregulated color: #0047AB

## Set colors for annotations
q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = "E") # E= Civis
v_colors

Gene_pHeat_Map <- function(results, title, padj = 0.01, lfc.abs = 1, sample.n = 1:8, c.width = c.width, c.height = c.height, f.size.row=5,
                           f.size.col=10, f.size=10, cols=colorRampPalette(c("#0047AB", "white", "#D22B2B"))(30), 
                           meta.annotation=c("Conditions", "Replicate"),
                           anno.colors = list(Replicate = c("1" = "#000004FF", "2" = "#5E626EFF"), 
                                              Conditions = c("0hr" = "#00204DFF", "12hr" = "#00336FFF", "48hr" = "#39486BFF", "144hr" = "#575C6DFF"))
) {
  require(pheatmap)
  require(viridis)
  require(dendsort)
  require(DESeq2)
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  ## Remove .## version of gene in ensembl ID
  #rownames(results) <- gsub("\\..*","",rownames(results))
  #rownames(dds) <- gsub("\\..*","",rownames(dds))
  
  
  ## Regularized log transformation, count extraction, and normalization
  rld <- DESeq2::rlog(dds, blind = FALSE)
  assay <-t(scale(t(assay(rld)[,sample.n])))
  
  ## Remove .## version of gene in ensembl ID
  rownames(assay) <- gsub("\\..*","",rownames(assay))
  
  ## ----
  res <- data.frame(res.48hr)
  sig.genes <- rownames(res[res$padj <= padj & abs(res$log2FoldChange) > abs(lfc.abs),])
  res <- as.data.frame(res[rownames(results) %in% sig.genes,])
  res <- res %>% 
    dplyr::distinct(symbol, .keep_all = T)
  
  res <- res[order(res$padj, decreasing = F),]
  
  res <- as.data.frame(res)[1:100,]
  res <-rownames(res)
  print(res)
  
  ## ----
  
  ## Organize genes by significance and extract
  sig.genes <- rownames(results[results$padj <= padj & abs(results$log2FoldChange) > abs(lfc.abs),])
  results <- as.data.frame(results[rownames(results) %in% sig.genes,])
  results <- results %>% 
    dplyr::distinct(symbol, .keep_all = T)
  
  ## Order genes by  pAdj, and extract N genes for heatmap
  results <- results[order(results$padj, decreasing = F),]
  #results <- as.data.frame(results[order(results$log2FoldChange, decreasing = FALSE),])[gene.n,]
  results <- as.data.frame(results[rownames(results) %in% res,])
  print(results)
  
  ## Prepare data frame for heatmapping 
  top.heat <- as.data.frame(assay[rownames(assay)%in%rownames(results),]) 
  top.heat
  ## Compute a distance calculation on both dimensions of the matrix for clustering
  #distance.gene <- dist(as.data.frame(assay[rownames(results),])) ## average distance method
  distance.gene <- as.dist((1 - cor(t(as.data.frame(assay[rownames(results),]))))/2) ## correlation method
  distance.sample <- dist(t(as.data.frame(assay[rownames(results),])))
  
  # Cluster based on the distance calculations
  #cluster.gene <- hclust(distance.gene, method="average") ## average distance method
  cluster.gene=hclust(distance.gene) ## correlation method
  cluster.sample <- hclust(distance.sample, method="average")
  
  ## Annotation information in data frame object
  dfcol <- as.data.frame(colData(dds)[,meta.annotation])
  
  # plot heat map
  p<-pheatmap(top.heat, 
              color = cols,
              fontsize = f.size,
              fontsize_row = f.size.row, 
              fontsize_col = f.size.col,
              cellwidth = c.width,
              cellheight = c.height,
              cluster_cols= F,
              #cluster_cols= cluster.sample,
              cluster_rows= cluster.gene,
              main = title,
              annotation_colors = anno.colors,
              labels_row = results$symbol,
              show_rownames = T,
              border_color = NA,
              annotation_legend = TRUE,
              annotation_col = dfcol
  )
  
  return(p)
  
}

## Use these for new plots
Gene_pHeat_Map(res.48hr.ov, "LAD overlapping", c.width = 9, c.height = 6) #48hr_inLADS_heatmap_top100
Gene_pHeat_Map(res.48hr.non, "LAD non-overlapping", c.width = 9, c.height = 6) #48hr_outsideLADS_heatmap_top100
##

#####################################################################################
# Gini Coefficients                                        
#####################################################################################

# upregulated color: #D22B2B
# downregulated color: #0047AB

library("DescTools")
library("ineq")
library("dplyr")
library(tidyr)
library(ggridges)
library(forcats)

##------------------------------Plot Histogram of Expression------------------------------##

## Break these plots up into deciles and see what it looks like 

plottheme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                 panel.grid.major = element_line(),
                 #panel.grid.minor = element_line(colour = "tomato", size=.25, linetype = "dashed"),
                 strip.background=element_blank(), axis.text.x=element_text(family="Arial", colour="dark blue"), axis.title.x=element_text(face="bold", size=20,family="Arial", colour="dark blue", vjust=-2),
                 axis.text.y=element_text(family="Arial", colour="dark blue"), 
                 axis.title.y=element_text(face="bold", size=15, family="Arial", colour="dark blue", vjust=2),
                 text = element_text(size=10, family="Arial"), 
                 plot.title=element_text(size=20, face="bold", family="Arial", color="tomato", hjust=0.0, vjust=10, lineheight=1.5),
                 plot.subtitle=element_text(size=15, face="bold", family="Arial", color="black", hjust=0.0, lineheight=1.5),
                 legend.title = element_text(size=12, color= "tomato",face="bold"), 
                 legend.text = element_text(size=10),
                 legend.key=element_rect(fill=NA),
                 axis.ticks=element_line(colour="black"))

## Remove .## version of gene in ensembl ID
rownames(dds) <- gsub("\\..*","",rownames(dds))

## Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Pre-filtering the dataset
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

## Get normalized counts
norm.counts <- counts(dds, normalized=T)
#norm.counts<- counts(dds)

## Scale data for plotting
#norm.counts <-scale(norm.counts) ## For gini calculations, don't use this, use unscaled counts
head(norm.counts)

## extract counts. Counts must be raw positive counts from non-normalized data
cts.wt <- norm.counts[,10:12] # 0hr counts (control)
cts.pol1 <- norm.counts[,4:6] # 12hr counts 
cts.actd <- norm.counts[,1:3] # 48hr counts 


## Get means of each group across replicates
mean.wt <- as.data.frame(apply(cts.wt , 1, function (x) mean(x)))
mean.pol1 <-as.data.frame(apply(cts.pol1, 1, function (x) mean(x)))
mean.actd <-  as.data.frame(apply(cts.actd, 1, function (x) mean(x)))

  
## Prepare data for plotting
mean.data <- cbind(mean.wt, mean.pol1, mean.actd)
colnames(mean.data) <- c("WT", "Pol 1", "ActD")

## Factor levels
factor1 <- rep(1, times = 1, length.out = nrow(mean.wt), each = 1)
factor2 <-rep(2, times = 1, length.out = nrow(mean.pol1), each = 1)
factor3 <-rep(3, times = 1, length.out = nrow(mean.actd), each = 1)
factor <- cbind(factor1, factor2, factor3)
factor <- as.data.frame(factor) %>% gather(Factor, value, na.rm = T)

## Rename with gene names and gather data into long format
mean.data$Genes <- rownames(norm.counts)
rownames(mean.data) <- make.names(rownames(norm.counts), unique=TRUE)
mean.data <- as.data.frame(mean.data) %>% gather(cond, Expression, -Genes, na.rm = T)
mean.data$Expression <- round(mean.data$Expression, 4)
mean.data$factor <- factor$value

## Set colors for individual plots
v_colors =  viridis(5, option = "E") # E= Civis
v_colors
max(mean.data$Expression)
mean(mean.data$Expression)
min(mean.data$Expression)

## Plot
p <- mean.data %>% 
  mutate(cond = fct_reorder(cond, factor)) %>% 
  ggplot(aes(x= cond, y= log10(Expression), fill = cond)) + 
  geom_violin(trim=F, size = 1) + theme_classic()+
  #scale_y_continuous(limits = c(0, 3))+ #c(3, 130)) for higher expressed genes
  #stat_summary(fun.data=data_summary, geom="pointrange", color="tomato")+
  scale_fill_manual(values=c("#414D6BFF", "#7C7B78FF", "#BCAF6FFF")) +
  labs(title="Scaled Expression", subtitle="", #subtitle="Top 1000 Variable Regions", 
       caption="")
p

mean.data <- mutate(mean.data, log10=log10(Expression))

ggdensity(mean.data, x = "log10",
          add = "mean", rug = TRUE,
          color = "cond", fill = "cond",
          palette = c("#414D6BFF", "#7C7B78FF", "#BCAF6FFF"))

##------------------------------Gini Coefficients------------------------------##

## Remember to use the unscaled counts, for all positive count values
## Function to filter count data for Gini calculation
filter_data <- function(data, filter){ #filter = lamin bed data, data = cts
  
  require(dplyr)
  
  data <- as.data.frame(data)
  data <-data %>%
    filter(rownames(data) %in% filter$ensembl_gene_id) 
  data <- as.data.frame(apply(data, 1, function (x) mean(x)))
  
  return(data)
}



## Filter data for each group
cts.0hr.non <- filter_data(cts.0hr, non.lamins)
cts.12hr.non <- filter_data(cts.12hr, non.lamins)
cts.48hr.non <- filter_data(cts.48hr, non.lamins)
cts.144hr.non <- filter_data(cts.144hr, non.lamins)

cts.0hr.ov <- filter_data(cts.0hr, lamin.overlaps)
cts.12hr.ov <- filter_data(cts.12hr, lamin.overlaps)
cts.48hr.ov <- filter_data(cts.48hr, lamin.overlaps)
cts.144hr.ov <- filter_data(cts.144hr, lamin.overlaps)

cts.0hr <- filter_data(cts.0hr, gene.pos)
cts.12hr <- filter_data(cts.12hr, gene.pos)
cts.48hr <- filter_data(cts.48hr, gene.pos)
cts.144hr <- filter_data(cts.144hr, gene.pos)

## Calculate Gini index for each treatment group
list = c(cts.0hr, cts.12hr, cts.48hr, cts.144hr, cts.0hr.ov,cts.12hr.ov, cts.48hr.ov, cts.144hr.ov, 
         cts.0hr.non, cts.12hr.non, cts.48hr.non, cts.144hr.non )

#list = c(cts.0hr, cts.48hr, cts.144hr, cts.0hr.ov,cts.48hr.ov, cts.144hr.ov, 
         #cts.0hr.non, cts.48hr.non, cts.144hr.non )

factor.list <- c()
gini.list <- c()
for (i in 1:length(list)) {
  this <- as.data.frame(list[i])
  this <- as.data.frame(apply(this, 2, function (x) ineq::Gini(x)))
  factor <- rep(i, times = 1, length.out = 1, each = 1)
  gini.list <-rbind(gini.list,this)
  factor.list <- c(factor.list, factor)
  
  
}

## Prepare data for plotting
rownames(gini.list) <- c("0hrs", "12hrs", "48hrs", "144hrs",
                         "0hrs in LADs", "12hrs in LADs", "48hrs in LADs", "144hrs in LADs",
                         "0hrs outside LADs", "12hrs outside LADs", "48hrs outside LADs", "144hrs outside LADs" )

#rownames(gini.list) <- c("0hrs", "48hrs", "144hrs",
                         #"0hrs in LADs", "48hrs in LADs", "144hrs in LADs",
                         #"0hrs outside LADs", "48hrs outside LADs", "144hrs outside LADs" )

gini.list <- cbind(gini.list, factor.list)
colnames(gini.list)  <- c("Gini", "factor")
gini.list$Time <- c("0hrs", "12hrs", "48hrs", "6 Day Washoff", "0hrs","12hrs", "48hrs", "6 Day Washoff", "0hrs","12hrs", "48hrs", "6 Day Washoff")
#gini.list$Time <- c("0hrs", "48hrs", "6 Day Washoff", "0hrs", "48hrs", "6 Day Washoff", "0hrs", "48hrs", "6 Day Washoff")

gini.list$factor <- c(1,2,3,4,1,2,3,4,1,2,3,4)
#gini.list$factor <- c(1,2,3,1,2,3,1,2,3)
gini.list$name <- c("All genes", "All genes", "All genes", "All genes", 
                    "Genes in LADs", "Genes in LADs", "Genes in LADs", "Genes in LADs",
                    "Genes outside LADs", "Genes outside LADs", "Genes outside LADs", "Genes outside LADs")
#gini.list$name <- c("All genes", "All genes", "All genes", 
                    #"Genes in LADs", "Genes in LADs", "Genes in LADs",
                    #"Genes outside LADs", "Genes outside LADs", "Genes outside LADs")

library(viridis)

## Set colors for individual plots
q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = "E") # E= Civis
v_colors


p <- gini.list %>% 
  mutate(Time = fct_reorder(Time, factor)) %>%
  ggplot(aes(x= Time, y= Gini, group=name)) + 
  geom_line(aes(color=name)) + geom_point(aes(color=name), size=4) + 
  #plottheme +
  scale_y_continuous( limits=c(0.8,1), breaks=c(0.8, 0.825,0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0))+
  scale_x_discrete(limits = c("0hrs","12hrs", "48hrs","6 Day Washoff"), expand = c(0.025, 0.025)) +
  scale_color_manual(values=c("#00204DFF", "#EE82EE", "#48D1CC"))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size = 13))+
  labs(title="Gini Coefficients", subtitle="Gini Coefficients", #subtitle="Top 1000 Variable Regions", 
       caption="")
p

# upregulated color: #D22B2B
# downregulated color: #0047AB

dev.off()

## Source code to calculate cumulutive sum (lorenz curve) taken from 'ineq' package
Lorenz <- function(x, n = rep(1, length(x)), plot = FALSE)
{
  ina <- !is.na(x)    
  n <- n[ina]
  x <- as.numeric(x)[ina]
  k <- length(x)
  o <- order(x)
  x <- x[o]
  n <- n[o]
  x <- n*x
  p <- cumsum(n)/sum(n) # if everything was distributed evenly
  L <- cumsum(x)/sum(x) # the distribution of the counts
  p <- c(0,p)
  L <- c(0,L)
  L2 <- L * mean(x)/mean(n)
  Lc <- list(p,L,L2)
  names(Lc) <- c("p", "L", "L.general")
  class(Lc) <- "Lc"
  if(plot) plot(Lc)
  Lc
}

## Get mean of each gene
treatment <-apply(cts.12hr, 2, function (x) mean(x))
control <-apply(cts.0hr, 2, function (x) mean(x))

## Apply lorenz function to gene averages for plotting
treat <-Lorenz(treatment)
treat$L
ctrl <-Lorenz(control)
ctrl$L

gene.color <- "#5302a3"
name <-"12 Hour"

## Plot in Base R
plot(treat, col = gene.color, main="",
     lty = 2, lwd = 2, xlab = "Replicate", 
     ylab = "Cumulative Expression") 
lines(ctrl, col = "#febb81", lty = 1, lwd = 2)

title(paste0("Lorenz curves, ", name, " , By Treatment"))

legend("topleft", lty = 1:2, lwd = 2, cex = 1.2, legend = 
         c("Perfect Distribution", name,
           "Control"),  
       col = c("black", gene.color, "#febb81")) 

#####################################################################################
# P10:P10                                 
#####################################################################################
nrow(mean.data)
mean.data <- mean.data[mean.data$Expression > 0.01, ]
nrow(mean.data)

p10.0hr <- mean.data[mean.data$Time == "0Hours",]
p10.0hr <- p10.0hr %>% mutate(Bin = ntile(Expression, n=10))
bottom.0hr <- p10.0hr[p10.0hr$Bin==1,]
top.0hr <- p10.0hr[p10.0hr$Bin==10,]
a <-mean(top.0hr$Expression)
b <- mean(bottom.0hr$Expression)
ratio1=abs(a/b)

p10.12hr <- mean.data[mean.data$Time == "12 Hours",]
p10.12hr <- p10.12hr %>% mutate(Bin = ntile(Expression, n=10))
bottom.12hr <- p10.12hr[p10.12hr$Bin==1,]
top.12hr <- p10.12hr[p10.12hr$Bin==10,]
a <-mean(top.12hr$Expression)
b <- mean(bottom.12hr$Expression)
ratio2=(a/b)

p10.48hr <- mean.data[mean.data$Time == "48 Hours",]
p10.48hr <- p10.48hr %>% mutate(Bin = ntile(Expression, n=10))
bottom.48hr <- p10.48hr[p10.48hr$Bin==1,]
top.48hr <- p10.48hr[p10.48hr$Bin==10,]
a <-mean(top.48hr$Expression)
b <- mean(bottom.48hr$Expression)
ratio3=abs(a/b)

p10.144hr <- mean.data[mean.data$Time == "144 Hours",]
p10.144hr <- p10.144hr %>% mutate(Bin = ntile(Expression, n=10))
bottom.144hr <- p10.144hr[p10.144hr$Bin==1,]
top.144hr <- p10.144hr[p10.144hr$Bin==10,]
a <-mean(top.144hr$Expression)
b <- mean(bottom.144hr$Expression)
ratio4=(a/b)

#####################################################################################
# LAD BED DEG intersection                                     
#####################################################################################

library(plyranges)
library(GenomicRanges)
library(tidyverse) 

#bed.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Raw_Data/BEDs"
bed.dir <- "/Volumes/My Passport/projects/Emily/Lamins_RNAseq/Raw_Data/BEDs"
list.files(bed.dir)

## Get chrom names and sizes (formatted for ensemble)
#chrom.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Raw_Data/Annotations"
chrom.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/101823_RNAseq/Annotations"
list.files(chrom.dir)

## Genomic coordinates for all coding genes in hg38
gene.pos <- read.table(file.path(chrom.dir, "hg38_gene_positions.csv"), sep= ",", header= T)

## Convert to Genomic Ranges data
gr.gene.pos <- as.data.frame(gene.pos)
gr.gene.pos <- gr.gene.pos %>% 
  transform( seqnames= paste0("chr",gr.gene.pos$chromosome_name), start = gr.gene.pos$start_position, end = gr.gene.pos$end_position)  %>% 
  as_granges()
gr.gene.pos

##-----------------------------------BED reps-----------------------------------##
## Using BED files from https://data.4dnucleome.org/experiments-damid/4DNEXMRJPX7X/ for Labmin B1
bed.rep1 <- read_bed(file.path(bed.dir, "4DNFI2BGIZ5F.bed") , col_names = NULL, genome_info = "hg38")
bed.rep1

bed.rep2 <- read_bed(file.path(bed.dir, "4DNFIC2T8L7T.bed") , col_names = NULL, genome_info = "hg38")
bed.rep2

bed.rep.union <- GenomicRanges::reduce(c(bed.rep1, bed.rep2))

##------------------------------Luay's Lamin B1 B2------------------------------##
## Using Luay's BED file for Lamin B1 and B2
bed.lmb1 <- read_bed(file.path(bed.dir, "LMB1_4DNFICCV71TZ.bed") , col_names = NULL, genome_info = "hg38")
bed.lmb1

bed.lmb2 <- read_bed(file.path(bed.dir, "LMB2_4DNFIBQH62LX.bed") , col_names = NULL, genome_info = "hg38")
bed.lmb2

## Union of lmb1 and lmb2 files
bed.lam.union <- GenomicRanges::reduce(c(bed.lmb1, bed.lmb2))
bed.lam.union

##----------------------------------Parse LADs----------------------------------##
## Name all of the peaks based on their genomic position
names(bed.lam.union) = paste0(seqnames(bed.lam.union),':',start(bed.lam.union),'-',end(bed.lam.union))

## Assign names to temp variable
region.names <- names(bed.lam.union)

## Find intersection between LADs and all genes
## intersect_rng <- join_overlap_intersect(query, subject) https://bioconductor.org/packages/devel/bioc/vignettes/plyranges/inst/doc/an-introduction.html
lamin.overlaps <- join_overlap_intersect(bed.lam.union, gr.gene.pos)

## Find non-overlapping regions https://support.bioconductor.org/p/74077/ sp over() method gr1[!gr1 %over% gr2,]
non.lamins <- gr.gene.pos[!gr.gene.pos %over% lamin.overlaps,]

##--------------------------------------------------##
##--If analyzing Lamin 1 and 2 overlaps seperately--##
lamin1.overlaps <- join_overlap_intersect(bed.lmb1, gr.gene.pos)
lamin2.overlaps <- join_overlap_intersect(bed.lmb2, gr.gene.pos)

##---------------------------Overlaps: Data Wrangling---------------------------##

## Convert to data frame
lamin.overlaps <- as.data.frame(lamin.overlaps, row.names = seq(1:29268))
##--------------------------------------------------##
##--If analyzing Lamin 1 and 2 overlaps seperately--##
#lamin1.overlaps <- as.data.frame(lamin1.overlaps, row.names = seq(1:22140))
#lamin2.overlaps <- as.data.frame(lamin1.overlaps, row.names = seq(1:28045))

## Remove non-unique Ensembl IDs
lamin.overlaps <-lamin.overlaps %>%  
  distinct(ensembl_gene_id,.keep_all = T ) 
rownames(lamin.overlaps) <-lamin.overlaps$ensembl_gene_id

## Convert results  to DF
res.12hr.df <- as.data.frame(res.12hr) 
res.48hr.df <- as.data.frame(res.48hr)
res.144hr.df <- as.data.frame(res.144hr) 

## Filter by ensembl ID and lfc
res.12hr.ov <-res.12hr.df %>%
  filter(rownames(res.12hr.df) %in% lamin.overlaps$ensembl_gene_id) 
res.48hr.ov <-res.48hr.df %>%
  filter(rownames(res.48hr.df) %in% lamin.overlaps$ensembl_gene_id) 
res.144hr.ov <-res.144hr.df %>%
  filter(rownames(res.144hr.df) %in% lamin.overlaps$ensembl_gene_id) 

## Volcano plots of lamin intersctions
DE_Vol_Plot(res.12hr.ov,"Lamin KO at 12 Hours: Within LADs", "res.12hr.LADs.volplot" )
DE_Vol_Plot(res.48hr.ov,"Lamin KO at 48 Hours: Within LADs", "res.48hr.LADs.volplot" )
DE_Vol_Plot(res.144hr.ov,"Lamin KO at 144 Hours: Within LADs", "res.144hr.LADsvolplot" )

## Volcano plots of lamin intersctions with p value =0.1 and lfc = 0.58
DE_Vol_Plot(res.12hr.ov,"Lamin KO at 12 Hours: Within LADs pv=0.1.", "res.12hr.LADs.pv0.1.volplot" )
DE_Vol_Plot(res.48hr.ov,"Lamin KO at 48 Hours: Within LADs pv=0.1.", "res.48hr.LADs.pv0.1.volplot" )
DE_Vol_Plot(res.144hr.ov,"Lamin KO at 144 Hours: Within LADs pv=0.1.", "res.144hr.LADs.pv0.1.volplot" )

##-------------------------Non Overlaps: Data Wrangling-------------------------##

## Convert to data frame
non.lamins <- as.data.frame(non.lamins, row.names = seq(1:33360))

## Remove non-unique Ensembl IDs
non.lamins <-non.lamins %>%  
  distinct(ensembl_gene_id,.keep_all = T ) 
rownames(non.lamins) <-non.lamins$ensembl_gene_id

## Filter by ensembl ID and lfc
res.12hr.non <-res.12hr.df %>%
  filter(rownames(res.12hr.df) %in% non.lamins$ensembl_gene_id) 
res.48hr.non <-res.48hr.df %>%
  filter(rownames(res.48hr.df) %in% non.lamins$ensembl_gene_id) 
res.144hr.non <-res.144hr.df %>%
  filter(rownames(res.144hr.df) %in% non.lamins$ensembl_gene_id)

## Volcano plots of non-lamin exluded genes
DE_Vol_Plot(res.12hr.non,"Lamin KO at 12 Hours: Outside LADs", "res.12hr.nonLADs.volplot" )
DE_Vol_Plot(res.48hr.non,"Lamin KO at 48 Hours: Outside LADs", "res.48hr.nonLADs.volplot" )
DE_Vol_Plot(res.144hr.non,"Lamin KO at 144 Hours: Outside LADs", "res.144hr.nonLADsvolplot" )

## Volcano plots of lamin intersctions with p value =0.1 and lfc = 0.58
DE_Vol_Plot(res.12hr.non,"Lamin KO at 12 Hours: Outside LADs pv=0.1.", "res.12hr.nonLADs.pv0.1.volplot" )
DE_Vol_Plot(res.48hr.non,"Lamin KO at 48 Hours: Outside LADs pv=0.1.", "res.48hr.nonLADs.pv0.1.volplot" )
DE_Vol_Plot(res.144hr.non,"Lamin KO at 144 Hours: Outside LADs pv=0.1.", "res.144hr.nonLADs.pv0.1.volplot" )

##-------------------------Box Plot: Lamin intersections-------------------------##

## Merge data
res.12hr.non <-res.12hr.non %>%
  mutate(group =rep("12hr nonoverlapping", nrow(res.12hr.non)))%>%
  mutate(factor =rep(6, nrow(res.12hr.non)))
res.12hr.ov <-res.12hr.ov %>%
  mutate(group =rep("12hr overlapping", nrow(res.12hr.ov)))%>%
  mutate(factor =rep(5, nrow(res.12hr.ov)))

res.48hr.non <-res.48hr.non %>%
  mutate(group =rep("48hr nonoverlapping", nrow(res.48hr.non)))%>%
  mutate(factor =rep(4, nrow(res.48hr.non)))
res.48hr.ov <-res.48hr.ov %>%
  mutate(group =rep("48hr overlapping", nrow(res.48hr.ov)))%>%
  mutate(factor =rep(3, nrow(res.48hr.ov)))

res.144hr.non <-res.144hr.non %>%
  mutate(group =rep("144hr nonoverlapping", nrow(res.144hr.non)))%>%
  mutate(factor =rep(2, nrow(res.144hr.non)))
res.144hr.ov <-res.144hr.ov %>%
  mutate(group =rep("144hr overlapping", nrow(res.144hr.ov)))%>%
  mutate(factor =rep(1, nrow(res.144hr.ov)))

## Plot all three conditions
plot.all <- rbind(res.12hr.non, res.12hr.ov, res.48hr.non, res.48hr.ov, res.144hr.non, res.144hr.ov)

## Plot just 12 and 144 hour conditions
plot.all <- rbind(res.12hr.non, res.12hr.ov, res.144hr.non, res.144hr.ov)

plot.all <-plot.all%>%
  filter(abs(log2FoldChange) > 1 ) %>%
  filter(padj < 0.01 )

plot.all <-plot.all %>%
  mutate(log2FoldChange=abs(log2FoldChange))

plot.all <-plot.all %>%
  mutate(padj=abs(log(padj)))

library(ggplot2)
library(forcats)

plottheme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                 panel.grid.major = element_line(),
                 panel.grid.minor = element_line(),
                 strip.background=element_blank(), axis.text.x=element_text(family="Arial", colour="dark blue"), axis.title.x=element_text(face="bold", size=20,family="Arial", colour="dark blue", vjust=-2),
                 axis.text.y=element_text(family="Arial", colour="dark blue"), 
                 axis.title.y=element_text(face="bold", size=15, family="Arial", colour="dark blue", vjust=2),
                 text = element_text(size=10, family="Arial"), 
                 plot.title=element_text(size=20, face="bold", family="Arial", color="tomato", hjust=0.0, vjust=10, lineheight=1.5),
                 plot.subtitle=element_text(size=15, face="bold", family="Aria", color="black", hjust=0.0, lineheight=1.5),
                 legend.title = element_text(size=12, color= "tomato",face="bold"), 
                 legend.text = element_text(size=10),
                 legend.key=element_rect(fill=NA),
                 axis.ticks=element_line(colour="black"))

## Box plot of absolute LFC for each group "#00204DFF", "#EE82EE", "#48D1CC"
plot.all%>% #Padj_Filter0.05_LADS_LFC_BoxPlot 12hrs_144hrs_Padj0.05_LADS_LFC_BoxPlot
mutate(group = fct_reorder(group, desc(factor))) %>%
ggplot(aes(x=as.factor(group), y=log2FoldChange, fill=group)) + 
  geom_boxplot(alpha=0.75, outlier.shape = NA) + plottheme + scale_fill_manual(values = c("#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE"))+
  xlab("") +  scale_y_continuous(name="Absolute LFC", limits=c(0.5,4), breaks=c(0.5,1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0))+
  labs(title="No P Adjusted Filter") + scale_x_discrete(labels=c("12Hrs outside LADs", "12Hrs inside LADs","48Hrs outside LADs", "48Hrs inside LADs","6 Days outside LADs", "6 Days inside LADs")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## Box plot of log P Adjusted for each group
plot.all %>% #Padj_Filter0.05_LADS_Padj_BoxPlot 12hrs_144hrs_Padj0.05_LADS_Padj_BoxPlot
  mutate(group = fct_reorder(group, desc(factor))) %>%
  ggplot(aes(x=as.factor(group), y=padj, fill=group)) + 
  geom_boxplot(alpha=0.75) + plottheme + scale_fill_manual(values = c("#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE"))+
  scale_y_continuous(name="Log Padj", limits=c(0,125), breaks=c(0,25, 50, 75, 100, 125))+
  labs(title="P Adjusted < 0.05") + scale_x_discrete(labels=c("12Hrs outside LADs", "12Hrs inside LADs","48Hrs outside LADs", "48Hrs inside LADs","6 Days outside LADs", "6 Days inside LADs")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


## Box plot of SEfor each group
plot.all %>% #Padj_Filter0.05_LADS_Padj_BoxPlot LADS_SE_BoxPlot
  mutate(group = fct_reorder(group, desc(factor))) %>%
  ggplot(aes(x=as.factor(group), y=lfcSE, fill=group)) + 
  geom_boxplot(alpha=0.75, outlier.shape = NA) + plottheme + scale_fill_manual(values = c("#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE"))+
  scale_y_continuous(name="Standard Error", limits=c(0,5), breaks=c(0,1,2,3,4,5))+
  labs(title="P Adjusted < 0.05") + scale_x_discrete(labels=c("12Hrs outside LADs", "12Hrs inside LADs","48Hrs outside LADs", "48Hrs inside LADs","6 Days outside LADs", "6 Days inside LADs")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


##------------------------Bar Plot: DEGs within/out LADs------------------------##

## Get DEGs for plotting relative fraction of DEGs in each group
res <- as.data.frame(res.144hr.non)
res <-res%>%
  filter(abs(log2FoldChange) > 1 )%>%
  filter(padj < 0.01 )

## Padj: 0.05 | LFC: 0.58
#dge <- c(2372,1813,551,4563,3471,1075,1902,1486,402)

## Padj: 0.01 | LFC: 1
dge <- c(559,144,413,1717,374,1331,412,104,296)
day <- c(rep("12 hours" , 3) , rep("48 hours" , 3) , rep("6 day washoff" , 3) )
cond <- rep(c("All DEG" ,"DEG in LADs", "DEG outside LADs") , 3)
factor <- c(9,8,7,6,5,4,3,2,1)
plot.df <- data.frame(cond, day, dge, factor)

## Plot box of DEGs
plot.df %>%
  mutate(day = fct_reorder(day, factor, .desc=T)) %>%
  ggplot(aes(y=dge, x=day, color = cond)) + plottheme+
  geom_bar(position=position_dodge(0.9), stat="identity", width=0.8, fill="white") + scale_color_manual(values = c("black", "#0047AB", "#D22B2B"))+
  geom_text(aes(label=dge), position = position_dodge(width = 0.875),hjust = 0.5, vjust=-0.5, size=3.5) + xlab("Day") + ylab("Differentially Expressed Genes")


# upregulated color: #D22B2B
# downregulated color: #0047AB

##----------------------Write Out LAD Intersection Results----------------------##

## Get DEGs for plotting relative fraction of DEGs in each group
res <- res.12hr.ov
res <-res%>%
  filter(abs(log2FoldChange) > 0.58 )%>%
  filter(padj < 0.05 )

dir = "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Results/"
filename = "12hrDEGs_WithinLADs"

write.csv(res, file=paste0(dir,filename, ".csv"), row.names=FALSE)

##----------------------Trajectory Plot----------------------##


## Organize genes by significance and extract
get_Ratio <- function(results.ov, results.non){
  degs.ov <- rownames(results.ov[results.ov$padj <= 0.01 & abs(results.ov$log2FoldChange) > abs(1),])
  degs.non <- rownames(results.non[results.non$padj <= 0.01 & abs(results.non$log2FoldChange) > abs(1),])
  ratio = as.numeric(round(length(degs.ov)/length(degs.non), 4))
  return(ratio)
}

hr.12 <- get_Ratio(res.12hr.ov, res.12hr.non)
hr.48 <- get_Ratio(res.48hr.ov, res.48hr.non)
hr.144 <- get_Ratio(res.144hr.ov, res.144hr.non)

hr.12.ov <- get_Ratio(res.12hr.ov, res.12hr)
hr.48.ov <- get_Ratio(res.48hr.ov, res.48hr)
hr.144.ov <- get_Ratio(res.144hr.ov, res.144hr)

hr.12.non <- get_Ratio(res.12hr.non, res.12hr)
hr.48.non <- get_Ratio(res.48hr.non, res.48hr)
hr.144.non <- get_Ratio(res.144hr.non, res.144hr)

list.non <- rbind(hr.12.non, hr.48.non, hr.144.non)
list.ov <- rbind(hr.12.ov, hr.48.ov, hr.144.ov)
list.all <- rbind(hr.12, hr.48, hr.144)

name <- c(rep(c("LAD DEGs/Non-LAD DEGs"),3), rep(c("LAD DEGs"),3),rep(c("Non-LAD DEGs"),3))
time <- rep(c("12hrs", "48hrs","6 Day Washoff"), 3)
factor <- rep(c(1,2,3), 3)
list <- rbind(list.all,list.ov, list.non)
list <- as.data.frame(cbind(list, name, time, factor))
colnames(list) <- c("ratio","name","time","factor")

#f3e79b,#fac484,#f8a07e,#eb7f86,#ce6693,#a059a0,#5c53a5

p <- list %>% 
  mutate(time = fct_reorder(time, factor)) %>%
  ggplot(aes(x= time, y= as.numeric(ratio), group = name)) + 
  geom_line(aes(color=name)) + geom_point(aes(color=name), size=8) + 
  plottheme +
  scale_y_continuous(limits = c(0,1),breaks=c(0.0, 0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8,0.9, 1.0))+
  scale_x_discrete(limits = c("12hrs", "48hrs","6 Day Washoff"), expand = c(0.025, 0.025)) +
 scale_color_manual(values=c("#f8a07e","#ce6693","#5c53a5"))+
 theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size = 15))+ ylab("DEG Ratio") + xlab("Timepoint") +
  ggtitle("Trajectory Plot")
p



#########################################################################################
# GO using Cluster Profiler
#########################################################################################


Enrich_GO <- function(results, ont, title, filename, n.results=10, x.axis="count", plot.col="qvalue"){
  
  require(clusterProfiler)
  require(org.Hs.eg.db)
  
  GO.res.dn <- as.data.frame(results[which(results$log2FoldChange <= -0.58 & results$padj <= 0.05),])
  GO.res.up <- as.data.frame(results[which(results$log2FoldChange >= 0.58 & results$padj <= 0.05),])

  egoup <- enrichGO( gene          = rownames(GO.res.up),
                     #universe      = names(results$symbol),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = ont,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable = TRUE)
  print(egoup)
  dot.p.up <- dotplot(egoup, x=x.axis, showCategory=n.results, color=plot.col ) + ggtitle(paste("UP",title))
  #ggsave(paste(dir.plt,"/",filename,"_UP",".png", sep=""), dpi=300)
  print(dot.p.up)
  
  egodn <- enrichGO(   gene          = rownames(GO.res.dn),
                       #universe      = names(results$symbol),
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = ont,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE)
  
  dot.p.dn <- dotplot(egodn, x=x.axis, showCategory=n.results, color=plot.col ) + ggtitle(paste("DN",title))
  print(dot.p.dn)
  #ggsave(paste(dir.plt,"/",filename,"_DN", ".png", sep=""), dpi=300)
}

Enrich_GO(res.pol1, "MF", " GO Enrichment MF: 12 Hours","12hr_vs_cDN_GO_MF_dotplot")
Enrich_GO(res.pol1, "CC", " GO Enrichment CC: 12 Hours","12hr_vs_cDN_GO_CC_dotplot")
Enrich_GO(res.pol1, "BP", " GO Enrichment BP: 12 Hours","12hr_vs_cDN_GO_BP_dotplot")

Enrich_GO(res.pol2, "MF", " GO Enrichment MF: 48 Hours","48hr_vs_cDN_GO_MF_dotplot")
Enrich_GO(res.pol2, "CC", " GO Enrichment CC: 48 Hours","48hr_vs_cDN_GO_CC_dotplot")
Enrich_GO(res.pol2, "BP", " GO Enrichment BP: 48 Hours","48hr_vs_cDN_GO_BP_dotplot")

Enrich_GO(res.actd, "MF", " GO Enrichment MF in LADs: 12 Hours","12hr_vs_cDN_GO_MF_inLADs_dotplot")
Enrich_GO(res.actd, "CC", " GO Enrichment CC in LADs: 12 Hours","12hr_vs_cDN_GO_CC_inLADs_dotplot")
Enrich_GO(res.actd, "BP", "GO Enrichment BP in LADs: 12 Hours","12hr_vs_cDN_GO_BP_inLADs_dotplot")


#####################################################################################
# LAD BED DEG intersection for LMBN1 and LMBN2                                    
#####################################################################################

library(plyranges)
library(GenomicRanges)
library(tidyverse) 

bed.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Raw_Data/BEDs"
list.files(bed.dir)

## Get chrom names and sizes (formatted for ensemble)
chrom.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Raw_Data/Annotations"
list.files(chrom.dir)

## Genomic coordinates for all coding genes in hg38
gene.pos <- read.table(file.path(chrom.dir, "hg38_gene_positions.csv"), sep= ",", header= T)

## Convert to Genomic Ranges data
gr.gene.pos <- as.data.frame(gene.pos)
gr.gene.pos <- gr.gene.pos %>% 
  transform( seqnames= paste0("chr",gr.gene.pos$chromosome_name), start = gr.gene.pos$start_position, end = gr.gene.pos$end_position)  %>% 
  as_granges()
gr.gene.pos

##------------------------------Luay's Lamin B1 B2------------------------------##
## Using Luay's BED file for Lamin B1 and B2
bed.lmb1 <- read_bed(file.path(bed.dir, "LMB1_4DNFICCV71TZ.bed") , col_names = NULL, genome_info = "hg38")
bed.lmb1

bed.lmb2 <- read_bed(file.path(bed.dir, "LMB2_4DNFIBQH62LX.bed") , col_names = NULL, genome_info = "hg38")
bed.lmb2

## Get unique regions in lmb1 and lmb2 #filter_by_non_overlaps(query, subject)
lmb1.unique <- bed.lmb1[!bed.lmb1 %over% bed.lmb2,]
lmb2.unique <- bed.lmb2[!bed.lmb2 %over% bed.lmb1,]

## Second way to get unique regions (plyranges)
lmb1.unique <- filter_by_non_overlaps(bed.lmb1, bed.lmb2)
lmb2.unique <- filter_by_non_overlaps(bed.lmb2, bed.lmb1)

##----------------------------------Parse LADs----------------------------------##
## Name all of the peaks based on their genomic position
names(lmb1.unique) = paste0(seqnames(lmb1.unique),':',start(lmb1.unique),'-',end(lmb1.unique))
names(lmb2.unique) = paste0(seqnames(lmb2.unique),':',start(lmb2.unique),'-',end(lmb2.unique))


## Assign names to temp variable
region.names.l1 <- names(lmb1.unique)
region.names.l2 <- names(lmb2.unique)

## Find intersection between LADs and all genes
## intersect_rng <- join_overlap_intersect(query, subject) https://bioconductor.org/packages/devel/bioc/vignettes/plyranges/inst/doc/an-introduction.html
lamin1.overlaps <- join_overlap_intersect(lmb1.unique, gr.gene.pos)
lamin2.overlaps <- join_overlap_intersect(lmb2.unique, gr.gene.pos)

##---------------------------Overlaps: Data Wrangling---------------------------##

## Convert to data frame
lamin1.overlaps <- as.data.frame(lamin1.overlaps, row.names = seq(1:627))
lamin2.overlaps <- as.data.frame(lamin2.overlaps, row.names = seq(1:439))

## Remove non-unique Ensembl IDs
lamin1.overlaps <-lamin1.overlaps %>%  
  distinct(ensembl_gene_id,.keep_all = T ) 
rownames(lamin1.overlaps) <-lamin1.overlaps$ensembl_gene_id

## Remove non-unique Ensembl IDs
lamin2.overlaps <-lamin2.overlaps %>%  
  distinct(ensembl_gene_id,.keep_all = T ) 
rownames(lamin2.overlaps) <-lamin2.overlaps$ensembl_gene_id

## Convert results  to DF
res.12hr.df <- as.data.frame(res.12hr) 
res.48hr.df <- as.data.frame(res.48hr)
res.144hr.df <- as.data.frame(res.144hr) 

## Filter by ensembl ID and lfc
res.12hr.l1 <-res.12hr.df %>%
  filter(rownames(res.12hr.df) %in% lamin1.overlaps$ensembl_gene_id) 
res.48hr.l1 <-res.48hr.df %>%
  filter(rownames(res.48hr.df) %in% lamin1.overlaps$ensembl_gene_id) 
res.144hr.l1 <-res.144hr.df %>%
  filter(rownames(res.144hr.df) %in% lamin1.overlaps$ensembl_gene_id) 

## Filter by ensembl ID and lfc
res.12hr.l2 <-res.12hr.df %>%
  filter(rownames(res.12hr.df) %in% lamin2.overlaps$ensembl_gene_id) 
res.48hr.l2 <-res.48hr.df %>%
  filter(rownames(res.48hr.df) %in% lamin2.overlaps$ensembl_gene_id) 
res.144hr.l2 <-res.144hr.df %>%
  filter(rownames(res.144hr.df) %in% lamin2.overlaps$ensembl_gene_id) 


##-------------------------Box Plot: Lamin intersections-------------------------##

## Merge data
res.12hr.l1 <-res.12hr.l1 %>%
  mutate(group =rep("12hr lamin 1", nrow(res.12hr.l1)))%>%
  mutate(factor =rep(6, nrow(res.12hr.l1)))
res.12hr.l2 <-res.12hr.l2 %>%
  mutate(group =rep("12hr lamin 2", nrow(res.12hr.l2)))%>%
  mutate(factor =rep(5, nrow(res.12hr.l2)))

res.48hr.l1 <-res.48hr.l1 %>%
  mutate(group =rep("48hr kamin 1", nrow(res.48hr.l1)))%>%
  mutate(factor =rep(4, nrow(res.48hr.l1)))
res.48hr.l2 <-res.48hr.l2 %>%
  mutate(group =rep("48hr lamin 2", nrow(res.48hr.l2)))%>%
  mutate(factor =rep(3, nrow(res.48hr.l2)))

res.144hr.l1 <-res.144hr.l1 %>%
  mutate(group =rep("144hr lamin 1", nrow(res.144hr.l1)))%>%
  mutate(factor =rep(2, nrow(res.144hr.l1)))
res.144hr.l2 <-res.144hr.l2 %>%
  mutate(group =rep("144hr lamin 2", nrow(res.144hr.l2)))%>%
  mutate(factor =rep(1, nrow(res.144hr.l2)))

## Plot all three conditions
plot.all <- rbind(res.12hr.l1, res.12hr.l2, res.48hr.l1, res.48hr.l2, res.144hr.l1, res.144hr.l2)

plot.all <-plot.all%>%
  filter(abs(log2FoldChange) > 0.58 ) %>%
  filter(padj < 0.1 )

plot.all <-plot.all %>%
  mutate(log2FoldChange=abs(log2FoldChange))

plot.all <-plot.all %>%
  mutate(padj=abs(log(padj)))

library(ggplot2)
library(forcats)

plottheme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                 panel.grid.major = element_line(),
                 panel.grid.minor = element_line(),
                 strip.background=element_blank(), axis.text.x=element_text(family="Arial",colour="dark blue"), axis.title.x=element_text(face="bold", size=20,family="Arial", colour="dark blue", vjust=-2),
                 axis.text.y=element_text(family="Arial",colour="dark blue"), 
                 axis.title.y=element_text(face="bold", size=15, family="Arial", colour="dark blue", vjust=2),
                 text = element_text(size=10, family="Arial"), 
                 plot.title=element_text(size=20, face="bold", family="Arial", color="tomato", hjust=0.0, vjust=10, lineheight=1.5),
                 plot.subtitle=element_text(size=15, face="bold", family="Aria", color="black", hjust=0.0, lineheight=1.5),
                 legend.title = element_text(size=12, color= "tomato",face="bold"), 
                 legend.text = element_text(size=10),
                 legend.key=element_rect(fill=NA),
                 axis.ticks=element_line(colour="black"))

## Box plot of absolute LFC for each group "#00204DFF", "#EE82EE", "#48D1CC"
plot.all%>% #Padj_Filter0.05_LADS_LFC_BoxPlot 12hrs_144hrs_Padj0.05_LADS_LFC_BoxPlot Padj0.1lamin1and2_LADS_LogPadj_BoxPlot
  mutate(group = fct_reorder(group, desc(factor))) %>%
  ggplot(aes(x=as.factor(group), y=log2FoldChange, fill=group)) + 
  geom_boxplot(alpha=0.75) + plottheme + scale_fill_manual(values = c("#48D1CC", "#EE82EE","#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE"))+
  xlab("") +  scale_y_continuous(name="Absolute LFC", limits=c(0.5,4), breaks=c(0.5,1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0))+
  labs(title="No P Adjusted Filter") + scale_x_discrete(labels=c("12Hrs Lamin1 LADs", "12Hrs Lamin2 LADs","48Hrs Lamin1 LADs", "48Hrs Lamin2 LADs","6 Days Lamin1 LADs", "6 Days Lamin2 LADs")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
## Box plot of log P Adjusted for each group
plot.all %>% #Padj_Filter0.05_LADS_Padj_BoxPlot 12hrs_144hrs_Padj0.05_LADS_Padj_BoxPlot
  mutate(group = fct_reorder(group, desc(factor))) %>%
  ggplot(aes(x=as.factor(group), y=padj, fill=group)) + 
  geom_boxplot(alpha=0.75) + plottheme + scale_fill_manual(values = c("#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE","#48D1CC", "#EE82EE"))+
  scale_y_continuous(name="Log Padj", limits=c(0,125), breaks=c(0,25, 50, 75, 100, 125))+
  labs(title="P Adjusted < 0.05") + scale_x_discrete(labels=c("12Hrs Lamin1 LADs", "12Hrs Lamin2 LADs","48Hrs Lamin1 LADs", "48Hrs Lamin2 LADs","6 Days Lamin1 LADs", "6 Days Lamin2 LADs")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


##------------------------Bar Plot: DEGs within/out LADs------------------------##

## Get DEGs for plotting relative fraction of DEGs in each group
res <- res.144hr.l2
res <-res%>%
  filter(abs(log2FoldChange) > 0.58 ) %>%
  filter(padj < 0.1 )

dge <- c(125,159,11,13,17,29,9,14)
day <- c( "All genes - Lamin1 LADs", "All genes - Lamin2 LADs", "12 hrs DEGs - Lamin1 LADs", "12 hrs DEGs - Lamin2 LADs", 
          "48 hrs DEGs - Lamin1 LADs", "48 hrs DEGs - Lamin2 LADs", "144 hrs DEGs - Lamin1 LADs", "144 hrs DEGs - Lamin2 LADs")

cond <- rep(c("Lamin1 LADs" , "Lamin2 LADs") , 4)
factor <- c(8,7,6,5,4,3,2,1)
plot.df <- data.frame(cond, day, dge, factor)

## Plot box of DEGs #DEGs_LADs_L1_L2
plot.df %>%
  mutate(day = fct_reorder(day, factor, .desc=T)) %>%
  ggplot(aes(y=dge, x=day, color = cond)) + plottheme+
  geom_bar(position=position_dodge(0.9), stat="identity", width=0.8, fill="white") + 
  scale_color_manual(values = c( "black","#0047AB", "#D22B2B"))+ theme(axis.text.x = element_text(angle = 45,  hjust=1) )+
  geom_text(aes(label=dge), position = position_dodge(width = 0.875),hjust = 0.5, vjust=-0.5, size=3.5) + xlab("Day") + ylab("Differentially Expressed Genes")

# upregulated color: #D22B2B
# downregulated color: #0047AB

##----------------------Write Out LAD Intersection Results----------------------##

## Get DEGs for plotting relative fraction of DEGs in each group
res <- res.144hr.l2
res <-res%>%
  filter(abs(log2FoldChange) > 0.58 )%>%
  filter(padj < 0.05 )

dir = "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Results/"
filename = "144hrDEGs_Within_LMB2_LADs"

write.csv(res, file=paste0(dir,filename, ".csv"), row.names=FALSE)

#####################################################################################
# Plot metascape results                                  
#####################################################################################

library(viridis)
library(forcats)

#OFvOP_RNAseq_Up_Stringent

dir_meta <-"/Volumes/My Passport/projects/Emily/Lamins_RNAseq/Results/Metascape"
list.files(dir_meta)
setwd(dir_meta)

## Read in Metascape results
meta.lamin <- read.table(file.path(dir_meta, "Metascape_Results_Summary.csv"), sep= ",", header= TRUE, fill = TRUE ) #OFvOP_RNAseq_Dn_Stringent

##  Build new data frame
meta.lamin <- data.frame(meta.lamin)

LAD = "in"
#LAD = "out"
Type = "cancer"
#Type = "structural"

meta.lamin <- meta.lamin[meta.lamin$LAD == LAD & meta.lamin$Type == Type,]

# upregulated color: #D22B2B
# downregulated color: #0047AB

subset = "12 Hours within LADs "
#subset = "12 Hours outside LADs "

## Plot categories # DisGeNET_within_Cancer
p.meta<-meta.lamin %>% 
  mutate(Name = fct_reorder(Name, pvalue)) %>% 
  ggplot(aes(x = pvalue, y = Name))+
  geom_point(aes(colour = as.numeric(pvalue),size = as.numeric(pvalue)))+
  scale_size_continuous(name="-Log 10 P Value")+
  scale_color_continuous(low = "#0047AB", high = "#D22B2B", name="P Value")+
  ggtitle(paste0("DisGeNET, ",subset," : ", Type))+ xlab("-Log 10 P Value")+ ylab("Annotation")+
  scale_x_continuous( limits = c(0,7), breaks = c(0,1,2,3,4,5,6,7))+
  theme(
    axis.line = element_line(color = "darkblue", 
                             size = 1, linetype = "solid"),
    axis.text.y = element_text( color="#008080", size=14, family = "Arial Narrow"),
    axis.text.x = element_text(angle = 0, color="black", size=12, family = "Arial Narrow"),
    plot.title = element_text(hjust = 0.5, color="black", face = "bold", size=15,
                              family = "Arial Rounded MT Bold"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank()) 
p.meta 




