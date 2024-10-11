#########################################################################################
## Gviz visualization ##
# used for Figure 7D-G & Figure 6H
# WARNING: code is very messy
#########################################################################################  
rm(list = ls())
gc()

library(biomaRt)
library(Gviz)
library(plyranges)
library(rtracklayer)
library(dplyr)
library(rtracklayer) 

options(ucscChromosomeNames=FALSE) ## turn off UCSC names

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

## BED files introduce these as annotation tracks as bed files
ctcf <- read.delim('/home/lmc0633/Transcription Paper/bigwigs/CTCF/ENCFF463FGL.bed',h = F, sep = '\t', stringsAsFactors = F)
rad21 <- read.delim('/home/lmc0633/Transcription Paper/bigwigs/Rad21/ENCFF391AAM.bed',h = F, sep = '\t', stringsAsFactors = F)
pol2ps5 <- read.delim('/home/lmc0633/Transcription Paper/bigwigs/Pol2ps5/ENCFF229YMU.bed',h = F, sep = '\t', stringsAsFactors = F)

pol2ps5<- pol2ps5%>% transform( seqnames= paste0(pol2ps5$V1), start = pol2ps5$V2, end = pol2ps5$V3) %>% as_granges()
rad21<- rad21%>% transform( seqnames= paste0(rad21$V1), start = rad21$V2, end = rad21$V3) %>% as_granges()
ctcf <- ctcf %>% transform( seqnames= paste0(ctcf $V1), start = ctcf $V2, end = ctcf $V3) %>% as_granges()

##-------------------------------------------------------## Bring in Pol 1 HMEC Chip

pol1 <- read.delim('/home/lmc0633/Transcription Paper/bigwigs/Pol1/GSM1544525_Pol_I_ChIP-seq_HMEC.txt',h = T, sep = '\t', stringsAsFactors = F)
colnames(pol1) <- c("seqnames", "start", "end")
pol1 <- na.omit(pol1)

a = 173607145-410000
b = 173616659+410000

y <- pol1 %>% filter(seqnames == 5)
gr<- y %>% 
  transform( seqnames= paste0(y$seqnames), start = y$start, end = y$end)  %>% 
  as_granges()

df <- data.frame(chrom = 5, start = a, end = b)
gr.2<- df %>% 
  transform( seqnames= chrom, start = start, end = end)  %>% 
  as_granges()

pol1 <- data.frame(gr %>% filter_by_overlaps(gr.2))

##-------------------------------------------------------## 

## BIGWIGS
k9me3.bw <- import.bw('/home/lmc0633/Transcription Paper/bigwigs/K9/ENCFF572IBD.bigWig',as="GRanges")
k4me3.bw <- import.bw('/home/lmc0633/Transcription Paper/bigwigs/K4/ENCFF394IMO.bigWig',as="GRanges") 

## RNAseq
actd.rna <- import.bw("/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Backman19_10.18.2023/totalRNA_BAMs/actd.merge.rpgc.bigWig",as="GRanges") 
wt.rna <- import.bw("/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Backman19_10.18.2023/totalRNA_BAMs/wt.merge.rpgc.bigWig",as="GRanges") 
pol1.rna <- import.bw("/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Backman19_10.18.2023/totalRNA_BAMs/pol1.merge.rpgc.bigWig",as="GRanges") 
pol2.rna <- import.bw("/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Backman19_10.18.2023/totalRNA_BAMs/pol2.merge.rpgc.bigWig",as="GRanges") 

## Proseq Pol2 AID2 data
ps.ctrl <- import.bw("/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Backman19_10.18.2023/totalRNA_BAMs/pol2_proseq/control/control.merged.rpgc.bigWig",as="GRanges") 
ps.ko <- import.bw("/projects/b1042/BackmanLab/Lucas/101823_rnaseq/Backman19_10.18.2023/totalRNA_BAMs/pol2_proseq/treatments/treatment.merged.rpgc.bigWig",as="GRanges") 

## transform scores to log2
actd.rna <- mutate(actd.rna, score= log2(score+1))
wt.rna <- mutate(wt.rna, score = log2(score+1))
pol1.rna <- mutate(pol1.rna, score = log2(score+1))
pol2.rna <- mutate(pol2.rna, score = log2(score+1))


seqlevelsStyle(ctcf) <- "NCBI"
seqlevelsStyle(rad21) <- "NCBI"
seqlevelsStyle(pol2ps5) <- "NCBI"
seqlevelsStyle(k4me3.bw) <- "NCBI"
seqlevelsStyle(k9me3.bw) <- "NCBI"
seqlevelsStyle(ps.ctrl) <- "NCBI"
seqlevelsStyle(ps.ko) <- "NCBI"
seqlevelsStyle(actd.rna) <- "NCBI"
seqlevelsStyle(wt.rna) <- "NCBI"
seqlevelsStyle(pol1.rna) <- "NCBI"
seqlevelsStyle(pol2.rna) <- "NCBI"

gc()

## Retrieve genes
values=c("ENSG00000147676") ## Mal2
values=c("ENSG00000145919") ## Bod1
values=c("ENSG00000177707") ## NECTIN3
values=c("ENSG00000023445") ## BIRC2 
values=c("ENSG00000286149") ## a set of TECs
values=c("ENSG00000169194") ## IL13

#gb <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','chromosome_name', "start_position","end_position","exon_chrom_start","exon_chrom_end", "strand"), filters =
#"ensembl_gene_id", values=values, mart=ensembl)

#colnames(gb) <- c('ENSG','ENST','chrom', 'txStart' , 'txEnd','exonStart' , 'exonEnd', 'strand')

#gb$chrom = gsub(gb$chrom, pattern = '^', replacement = 'chr')
#gb$strand = ifelse(gb$strand == 1, '+',"-")

## Plot gene track 
biomTrack <- BiomartGeneRegionTrack(genome = "hg38",gene=values, name = "ENSEMBL", biomart = ensembl, fill = "black", color="black", cex.feature = 0.75,
                                    fontsize=12,fontcolor.title="white",background.title = "black", fontcolor.group="black")

##-- Plotting bigwig tracks --##

chrom = 11
## Colors
# green "#798234"
# K9 red "#cf597e"
# K4 blue "#3969AC"
# K27 orange "#E68310"
#"#f3e79b","#fac484","#f8a07e","#eb7f86","#ce6693","#a059a0","#5c53a5"

##

dTrack1 <- DataTrack(range = k9me3.bw, genome = "hg38",
                     chromosome = chrom, 
                     name = "H3K9me3", type = "hist", 
                     window = -1, windowSize = 500, 
                     ylim=c(1,5),col.histogram="#ce6693", background.title = "#ce6693")

dTrack2 <- DataTrack(range = k4me3.bw, genome = "hg38",  
                     chromosome = chrom, 
                     name = "H3K4me3", type = "hist",
                     window = -1, windowSize = 100, 
                     ylim=c(1,30), col.histogram="#5c53a5", background.title = "#5c53a5")

##-------------------------------------------------------## Proseq tracks
ps.ctrl.track <- DataTrack(range = ps.ctrl, genome = "hg38",
                           chromosome = chrom, 
                           name = "PS CTRL", type = "hist", 
                           window = -1, windowSize = 1000, 
                           ylim=c(0,10),col.histogram="black", background.title = "black")

ps.ko.track  <- DataTrack(range = ps.ko, genome = "hg38",  
                          chromosome = chrom, 
                          name = "PS KO", type = "hist",
                          window = -1, windowSize = 1000, 
                          ylim=c(0,10), col.histogram="black", background.title = "black")

##-------------------------------------------------------## Proseq tracks

#aTrack3 <- AnnotationTrack(range = rad21, genome = "hg38", background.title = "black",name = "RAD21", chromosome = chrom, fill="#f3e79b", col ="#f3e79b" )

#aTrack4 <- AnnotationTrack(range = ctcf, genome = "hg38", background.title = "black",name = "CTCF", chromosome = chrom, fill="#f3e79b", col ="#f3e79b")

aTrackpol1 <- AnnotationTrack(range = pol1, genome = "hg38", background.title = "black",name = "POL1", chromosome = chrom, fill="black", col ="black")

aTrack5 <- AnnotationTrack(range = pol2ps5, genome = "hg38", background.title = "black",name = "POL2-PS5", chromosome = chrom, fill="#f8a07e", col ="#f8a07e")

dTrack6 <- DataTrack(range = pol1.rna, genome = "hg38",
                     chromosome = chrom, 
                     name = "Pol 1", type = "hist", 
                     window = -1, windowSize = 1000,
                     ylim=c(0,5),col.histogram="#a059a0", background.title = "black")

dTrack7 <- DataTrack(range = pol2.rna,type = "hist",genome = "hg38",
                     chromosome = chrom, 
                     window = -1, windowSize = 1000, ylim=c(0,5),
                     col.histogram="black", background.title = "black",
                     name = "Pol 2")

dTrack8 <- DataTrack(range = actd.rna, genome = "hg38",
                     chromosome = chrom,
                     name = "ActD", type = "hist", 
                     window = -1, windowSize = 1000,
                     ylim=c(0,5),col.histogram="black", background.title = "black")

dTrack9 <- DataTrack(range = wt.rna, genome = "hg38",name = "WT",
                     chromosome = chrom, type = "hist",window = -1, windowSize = 1000,
                     ylim=c(0,5), col.histogram="black", background.title = "black")



## MAL2 8 119165034 119165034
## BOD1 5	173607145	173616659
## NECTIN3 3	111070071	111275563 ## track.nectin3
##	BIRC3		11	102317484	102339403
## TECs 21	8462320	8462798
## IL13/Rad50 5 132526809-132678639
x = 102317484-150000
y = 102339403+150000

x=132526809-25000
y=132678639+25000

## Plot full locus
plotTracks(list( biomTrack,dTrack1, dTrack2, aTrackpol1, aTrack5, dTrack6, dTrack7, dTrack8, dTrack9), from = x, to =y,  collapseTranscripts="meta", showID=T, transcriptAnnotation="symbol")

## Plot Proseq
plotTracks(list( biomTrack, ps.ko.track, ps.ctrl.track), from = x, to =y,  collapseTranscripts="meta", showID=T,collapseTranscripts="meta", transcriptAnnotation="symbol")

#---------------------------------------------------------------#
# Author: Lucas Carter                                                   
# Email: lucascarter2025@u.northwestern.edu                                                          
# PI: Vadim Backman                                                   
# Description: 
# This script is for generating BigWig tracks using
# our Cut&Tag data
#
##-------------------------------------------------------## load packages

rm(list = ls())
gc()

library(biomaRt)
library(Gviz)
library(plyranges)
library(dplyr)
require(tidyr)
require(ggplot2)
require(viridis)
require(ggpubr)
library(rtracklayer)

##-------------------------------------------------------##  rDNA tracks plotting

## load BIGWIGS
bw.path <- '/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/mergedR/coverage'
list.files(bw.path)

p1.bw <- file.path(bw.path, "P1AB_WT.cpm.bigWig")
p2.bw <- file.path(bw.path, "P2AB_WT.cpm.bigWig")
k4.bw <- file.path(bw.path, "K4AB_WT.cpm.bigWig")
k9.bw <- file.path(bw.path, "K9AB_WT.cpm.bigWig")

## path to annotation BED
ann.path <- '/projects/b1042/BackmanLab/Lucas/090124_CnT/genome/refR/Human_hg38-rDNA_genome_v1.0_annotation'
ann.bed <- import.bed(file.path(ann.path, 'hg38-rDNA_v1.0.bed'))
ann.bed <- ann.bed %>% filter(seqnames =="chrR") ## restrict ranges

## plot BIGWIG tracks 
dTrack1 <- DataTrack(range = p1.bw, genome = "hg38",
                     chromosome = "chrR", 
                     name = "Pol 1", type = "histogram", 
                     window = -1, windowSize = 10, 
                     ylim=c(0,125),col.histogram="#fde725", col.sampleNames = "white", background.title = "black", fontcolor.group="white")

dTrack2 <- DataTrack(range = p2.bw, genome = "hg38",
                     chromosome = "chrR", 
                     name = "Pol 2", type = "histogram",
                     window = -1, windowSize = 10, 
                     ylim=c(0,200), col.histogram="#5ec962", background.title = "black")

dTrack3 <- DataTrack(range = k4.bw, genome = "hg38",
                     chromosome = "chrR", 
                     name = "K4Me3", type = "histogram", 
                     window = -1, windowSize = 10, 
                     ylim=c(0,5),col.histogram="#21918c", background.title = "black")

dTrack4 <- DataTrack(range = k9.bw, genome = "hg38",
                     chromosome = "chrR", 
                     name = "K9Me3", type = "histogram", 
                     window = -1, windowSize = 500,
                     ylim=c(0,70),col.histogram="#3b528b",fontcolor.title="black",background.title = "white", color="black")



## Filter out extra annotations
ann.ls <- c("IGS", "Spacer_Promoter", "Enhancer_Repeats", "47S_Promoter", "5'_ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'_ETS", "IGS")
ann.ls <- c("IGS", "18S", "5.8S", "28S", "IGS")

ann.bed <- ann.bed %>% filter(ann.bed$name %in% ann.ls)

## add in 
txTr <- GeneRegionTrack(ann.bed, chromosome = "chrR", start = 1,  end = 40000, symbol = ann.bed$name, name = "rDNA",fill = "black", color="black", cex.feature = 0.75,
                        fontsize=12,fontcolor.title="white",background.title = "black", fontcolor.group="black")

## Plot all tracks 
plotTracks(list(txTr, dTrack1, dTrack2, dTrack3, dTrack4),from = 1, to = 40000, showID=T,transcriptAnnotation="symbol")

