#########################################################################################
## HIC Compare ##
# This code contains everything necessary to regenerate the
# multiHiCcompare differential interaction analysis results
# DI analysis is very memory heavy and like other HiC data analysis
# must be done on HPC
#########################################################################################

rm(list = ls())
gc()

library(multiHiCcompare)
library(BiocParallel)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)

## set paralell processing
numCores <- 50
numCores <- parallel::detectCores()
register(MulticoreParam(workers = numCores-1), default = TRUE)

##-------------------------------------------------------## parse contacts for MultiHiCcompare  first

## loop for import
getContacts <- function(exps, rep, dir.contacts = "/projects/p32171/HiC2/opt/juicer/work/112123_HiC/contact_data"){
  
  files <- list.files(paste0(dir.contacts,"/",exps[1],"/",rep,"/intracontacts/"), pattern=".csv", all.files=T, include.dirs = FALSE)
  files<- files
  
  df <- data.frame()
  for (exp in 1:length(exps)){
    
    exp <- exps[exp]
    
    filepath <- paste0(dir.contacts,"/",exp,"/",rep,"/intracontacts/")
    cat(paste0("input path: ",filepath, "\n","\n"))
    
    
    for (file in 1:length(files)){
      infile <- files[file]
      
      name =sapply(strsplit(infile, "_"), `[`, 1)
      chrom=sapply(strsplit(name, "-"), `[`, 2)
      chrom=sapply(strsplit(chrom, "chr"), `[`, 2)
      chrom<- ifelse(chrom== "X", as.numeric(23), as.numeric(chrom))
      
      cat(paste('chromosome ',chrom,' condition',exp, '\n\n'))
      
      cat(paste('appending contact files',exp,":",infile,'\n\n'))
      
      infile <- read.table(file.path(filepath, infile), sep= ",", header= T)
      infile <-data.frame(infile)
      infile <- na.omit(infile)
      
      ## Get chroms and conds
      infile$chr <- rep(paste0(chrom),nrow(infile))
      
      df <- rbind(df, infile)
      
    }
    
  }
  
  colnames(df) <- c("region1", "region2", "IF", "chr")
  df <- df[, c("chr", "region1", "region2", "IF" )]
  
}

exps <- c("1hr_ActD",  "Pol2_6hrs_Aux","Pol1_6hrs_Aux", "WT_HCT116_CTRL")

pol1.r1 <- getContacts(exps = exps[3], rep="Rep1")
pol1.r2 <- getContacts(exps = exps[3], rep="Rep2")
wt.r1 <- getContacts(exps = exps[4],rep="Rep1")
wt.r2 <- getContacts(exps = exps[4],rep="Rep2")
pol2.r1 <- getContacts(exps = exps[2], rep="Rep1")
pol2.r2 <- getContacts(exps = exps[2], rep="Rep2")

dir <- "/projects/p32171/HiC2/opt/juicer/work/112123_HiC/contact_data/multiHiCompare/"
# Save an object to a file
saveRDS(pol1.r1, file = paste0(dir,"pol1_r1_contacts.rds"))
saveRDS(pol1.r2, file = paste0(dir,"pol1_r2_contacts.rds"))
saveRDS(pol2.r1, file = paste0(dir,"pol2_r1_contacts.rds"))
saveRDS(pol2.r2, file = paste0(dir,"pol2_r2_contacts.rds"))
saveRDS(wt.r1, file = paste0(dir,"wt_r1_contacts.rds"))
saveRDS(wt.r2, file = paste0(dir,"wt_r2_contacts.rds"))

## Put experiments into list
hics.list <- list(a = wt.r1,b = wt.r2,c = pol1.r1,d =pol1.r2)

##-------------------------------------------------------## MultiHiCcompare first route to analysis (loop to process all conditions)

a = pol1.r1
b = pol1.r2
c = pol2.r1
d = pol2.r2

exp= "Pol1vsPol2"
FDR_threshold <- 0.25
logFC_threshold <- 0.5

chroms <- paste0(seq(1,23,1))

for (chrom in 1:length(chroms)){
  
  chrom = chroms[chrom]
  
  
  w <- a[a$chr == chrom,]
  x <- b[b$chr== chrom,]
  y <- c[c$chr== chrom,]
  z <- d[d$chr== chrom,]
  
  hic.list <- list(w,x,y,z)
  cat(paste0("hic list generated ","\n"))
  
  cat(paste0("getting DIs for: ", chrom))
  
  DI <- hic.list |> 
    make_hicexp(
      data_list = hic.list, 
      groups = factor(c(0,0,1, 1)),  zero.p = 0.8, A.min = 5, filter = TRUE
    ) |> 
    fastlo(parallel = T) |> 
    hic_exactTest(p.method = 'fdr',parallel = T) |> 
    results()
  DI
  
  DI <- DI %>% filter(p.adj < FDR_threshold, abs(logFC) > logFC_threshold)
  
  dir.data <- "/projects/p32171/HiC2/opt/juicer/work/112123_HiC/contact_data/multiHiCompare/"
  file.name <- file.path(dir.data, paste0(exp,"_", chrom,"_topDI.txt"))
  write.table(data.frame(DI),file.name,sep="\t",col.names=F,row.names=F,quote=F)
  
}

gc()
##-------------------------------------------------------## MultiHiCcompare second route to analysis

#hicexp <- make_hicexp(wt.r1,wt.r2, pol1.r1, pol1.r2, groups = c(0,0, 1,1), zero.p = 0.8, A.min = 5, filter = TRUE)

## differential interaction with MultiHiCcompare with full cyclic loess (slow)          
DI <- hics.list |> 
  make_hicexp(
    data_list = hics.list, 
    groups = factor(c(0,0,1, 1)),  zero.p = 0.8, A.min = 5, filter = TRUE
  ) |> 
  cyclic_loess(parallel = T) |> 
  hic_exactTest(p.method = 'fdr',parallel = T) |> 
  results()
DI

## differential interaction with MultiHiCcompare with partial cyclic loess (fast)     
DI <- hics.list |> 
  make_hicexp(
    data_list = hics.list, 
    groups = factor(c(0,0,1, 1)),  zero.p = 0.8, A.min = 5, filter = TRUE
  ) |> 
  fastlo(parallel = T) |> 
  hic_exactTest(p.method = 'fdr',parallel = T) |> 
  results()
DI

FDR_threshold <- 0.25
logFC_threshold <- 0.5

DI <- DI %>% filter(p.adj < FDR_threshold, abs(logFC) > logFC_threshold)

exp= "Pol1"
dir.data <- "/home/lmc0633/Transcription Paper"
file.name <- file.path(dir.data, paste0(exp,"_topDI_TEST.bedpe"))
write.table(data.frame(DI),file.name,sep="\t",col.names=F,row.names=F,quote=F)

##-------------------------------------------------------## Write out results

# Path with data
dir.data <- "/home/lmc0633/Transcription Paper"

# File name to save the results
file.name <- file.path(dir.data, paste0(exp,"_DI_table.csv"))

# Save the data
write.csv(as.data.frame(DI), file.name)
list.files(dir.data)

## Bring in DI object
DI <- read.table(file.path(dir.data, "100kb_Reps_DI_table.csv"), sep= ",", header= T)
DI <- read.table(file.path(dir.data, "100kb_Pol2_DI_table.csv"), sep= ",", header= T)
DI <- read.table(file.path(dir.data, "100kb_actd_DI_table.csv"), sep= ",", header= T)

##-------------------------------------------------------## Volcano plot of DI results

## add a column of NAs
DI$sign <- "NO"
## If log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
DI$sign[DI$p.value <= 0.05 & DI$logFC >= 1] <- "UP"
## If log2Foldchange < -0.6 and pvalue < 0.05, set as "DN"
DI$sign[DI$p.value <= 0.05 & DI$logFC <= -1] <- "DN"

ggplot(DI, aes(x = logFC, y = -log10(p.value), col = sign)) + 
  geom_point(size = 0.3) + 
  theme_bw() + labs(subtitle=paste(paste("Downregulated: ",sum(DI$sign == "DN")), paste("| Upregulated: ",sum(DI$sign == "UP")))) +
  ylim(c(0, 40)) + geom_vline(xintercept=c(-1, 1), col="black", lty = "dashed") +
  geom_hline(yintercept=-log10(0.05) , col="black", lty = "dashed")+ 
  theme(legend.position = 'none') +   theme( legend.position="none",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(size = 8)) +
  #scale_color_manual(values = c("#D46780",'grey', "#798234"))
  scale_color_manual(values = c("black",'lightgrey', "black"))

##-------------------------------------------------------## MA plot of DI results

## Downsample data and r
set.seed(10)
plot <- DI[sample(nrow(DI), 2000000), ]

plot <- data.frame(baseMean = plot$D, log2FoldChange = plot$logFC, padj = plot$p.value)

library(ggpubr)

ggmaplot(plot, 
         fdr = 0.10, fc = 1, size = 2.5,
         #palette = c("#798234", "#D46780", "lightgray"),
         palette = c("black", "black", "lightgray"),
         legend = "top", top = 0,
         alpha=1,
         xlab = "Log10 Distance",
         ylab = "Log2 Fold Change",
         font.legend = c(16),
         ggtheme = theme_classic()
) +  theme( legend.position="top",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(size = 16)) +
  scale_y_continuous(limits = c(-7,7),breaks = seq(-7,7,1),labels = seq(-7,7,1)) +geom_hline(yintercept= c(-1, 1), col="black", lty = "dashed")

#########################################################################################
# Cis Trans ratios (generated using multiHiCcompare)
#########################################################################################

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(forcats)

list.files("/home/lmc0633/Transcription Paper")

ct.wt <- fread(file.path("/home/lmc0633/Transcription Paper", "ct.wt.csv"), sep= ",", header= T)
ct.pol1 <- fread(file.path("/home/lmc0633/Transcription Paper", "ct.pol1.csv"), sep= ",", header= T)
ct.actd <- fread(file.path("/home/lmc0633/Transcription Paper", "ct.actd.csv"), sep= ",", header= T)
ct.pol2 <- fread(file.path("/home/lmc0633/Transcription Paper", "ct.pol2.csv"), sep= ",", header= T)

## Get rid of other chroms
chroms <- paste0(c(as.character(1:22), "X"))
ct.wt <- ct.wt[ct.wt$chr %in% chroms,]
ct.pol1 <- ct.pol1[ct.pol1$chr %in% chroms,]
ct.pol2 <- ct.pol2[ct.pol2$chr %in% chroms,]
ct.actd <- ct.actd[ct.actd$chr %in% chroms,]

## Label conds
ct.wt$cond <- rep("WT", nrow(ct.wt))
ct.pol1$cond <- rep("POLR1A", nrow(ct.pol1))
ct.pol2$cond <- rep("POLR2A", nrow(ct.pol2))
ct.actd$cond <- rep("ActD", nrow(ct.actd))

## Reformat data for plotting
ct.plot <- rbind(ct.wt, ct.pol1)
ct.plot <- rbind(ct.wt, ct.pol2)
ct.plot <- rbind(ct.wt, ct.actd)
ct.plot <- rbind(ct.pol1, ct.pol2)
ct.plot<- ct.plot[,-c(1,3:5)]

ct.plot <- gather(ct.plot, key="label", value="pct", 2:3)
ct.plot <- mutate(ct.plot, chr = paste0(cond," chr",chr))

## Factor
ct.plot$factor <- rep(c(1:23),2)

## Plot 
ct.plot %>% 
  mutate(chr = fct_reorder(chr, as.numeric(factor))) %>% 
  ggplot( aes(fill=label, y=pct, x=chr, group = cond)) + 
  geom_bar(position="stack", stat="identity") +
  guides(x=guide_axis(angle = 90)) + theme_classic()+
  scale_y_continuous(labels = scales::percent) + coord_cartesian(ylim = c(0.5, 1))+
  theme( legend.position="none",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(size = 8)) +
  scale_fill_manual(values = c("darkgrey", "black"))+
  labs(x = 'Chromosomes', y = 'Cis Trans Ratio')

