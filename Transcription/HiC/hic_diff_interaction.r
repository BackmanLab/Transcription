#########################################################################################
# HIC Compare
#########################################################################################

rm(list = ls())
gc()

library(multiHiCcompare)
library(BiocParallel)
library(data.table)
library(tidyr)
library(dplyr)

## set paralell processing
numCores <- 16
register(MulticoreParam(workers = numCores-1), default = TRUE)
print(numCores)

## loop for import
getContacts <- function(exps, rep, dir.contacts = "/projects/p32171/HiC2/opt/juicer/work/112123_HiC/contact_data"){

files <- list.files(paste0(dir.contacts,"/",exps[1],"/",rep,"/intracontacts/"), pattern=".csv", all.files=T, include.dirs = FALSE)
files
  
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

pol2.r1 <- getContacts(exps = exps[2], rep="Rep1")
pol2.r2 <- getContacts(exps = exps[2], rep="Rep2")
pol1.r1 <- getContacts(exps = exps[3],rep="Rep1")
pol1.r2 <- getContacts(exps = exps[3],rep="Rep2")

## Put experiments into list
hics.list <- list(pol2.r1,pol2.r2,pol1.r1,pol1.r2)

## differential interaction with MultiHiCcompare            
DI <- hics.list |> 
  make_hicexp(
    data_list = hics.list, 
    groups = factor(c(1, 2)),  zero.p = 0.8, A.min = 5, filter = TRUE
  ) |> 
  cyclic_loess(parallel = T) |> 
  hic_exactTest(p.method = 'fdr',parallel = T) |> 
  results()
DI

##-------------------------------------------------------## export data
FDR_threshold <- 0.2
logFC_threshold <- 0.5
logCPM_threshold <- 0.5
td <- topDirs(DI, 
              logfc_cutoff = logFC_threshold, 
              logcpm_cutoff = logCPM_threshold,
              p.adj_cutoff = FDR_threshold, 
              return_df = 'pairedbed')

exp= "Pol1vsPol2"
dir.data <- "/home/lmc0633/Transcription Paper"
file.name <- file.path(dir.data, paste0(exp,"_topDI.bedpe"))
write.table(td,file.name,sep="\t",col.names=F,row.names=F,quote=F)
