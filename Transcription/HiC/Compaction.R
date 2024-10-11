#########################################################################################
## Compaction analysis by Loci ##
# This script requires contacts from .hic dumped using straw
# change directories to your dir where contacts are
#########################################################################################

library(data.table)
library(dplyr)
library(tidyr)

## Condition labels
exps <- c("WT_HCT116_CTRL", "Pol1_6hrs_Aux", "Pol2_6hrs_Aux",  "1hr_ActD")

# Set directory variables here
dir.contacts <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data" # our contact data
list.files(dir.contacts)

files <- list.files(paste0(dir.contacts,"/",exps[1],"/mega/intracontacts/"), pattern=".csv", all.files=T, include.dirs = FALSE)
files

##-------------------------------------------------------## Fine Grain Compaction Analysis Function

## Moving window contact log2 ratio for fine grain compaction analysis
mwa.cont.ratio <- function(distal.1, distal.2, local.1, local.2,
                           inc, window, 
                           dir, exps,
                           cols, col.num) {
  
  ## Initial empty data frame
  tbl<- data.frame()
  
  ## loop for calculating sliding window log2 ratio
  for (exp in 1:length(exps)){
    
    ## Loop through conditions
    exp <- exps[exp]
    files <- list.files(paste0(dir,"/",exp,"/mega/intracontacts/"), pattern=".csv", all.files=T, include.dirs = FALSE)
    
    for (file in 1:length(files)){
      
      ## Loop through contact files in each condition
      infile <- files[file]
      
      
      infile <- fread(file.path(paste0(dir.contacts,"/",exp,"/mega/intracontacts/"), infile), sep= ",", header= T)
      infile <- mutate(infile, width = abs(infile$x- infile$y)) ## Get distance of contact
      
      ## Grab chrom names
      name=sapply(strsplit(files[file], "_"), `[`, 1)
      chrom=sapply(strsplit(name, "-"), `[`, 2)
      chrom=sapply(strsplit(chrom, "chr"), `[`, 2)
      
      cat(paste('gathering log2 ratios for condition: ',exp," & Chromosome: ",chrom,'\n\n'))
      
      ## Generate list spaced by increment and length of infile left contacts
      list <- seq(min(infile$x),max(infile$x), by=inc) 
      
      for (i in (min(list)/inc):length(list)) {
        
        ## Collect rows between the size of increment and increment + window
        i = i*inc
        mwa <- infile[infile$x >= (i) & infile$x <= (i + window),]
        #print(paste0("number of rows: ", nrow(mwa)))
        
        wind <- paste0("window: ", chrom, " ",i, " - ",(i + window), " BP") ## Window string
        
        ## gather distal and local dataframes
        distal.int <- mwa[as.numeric(mwa$width)  >= distal.1 & as.numeric(mwa$width)  <= distal.2,] #distal.1 = lower bound int, distal.2 = upper bound int
        local.int <- mwa[as.numeric(mwa$width)  >= local.1 & as.numeric(mwa$width) <= local.2 ,] #local.1 = lower bound int, local.2 = upper bound int
        
        ## sum contact frequencies for distal and local groups 
        distal <-  sum(distal.int$counts,na.rm=T)
        local <-  sum(local.int$counts, na.rm=T)
        
        ## Calculate log ratio
        ratio = log2((distal/local))
        df <- data.frame(exp, chrom, i, (i + window), wind, ratio)
        
        tbl <- rbind(tbl,df)
        
      }
      
      
    }
    
    cat(paste(exp," Log2 DTL is done processing ",'\n\n'))
    
  }
  
  colnames(tbl) <- c("cond", "chrom", "start", "end", "window", "ratio")
  
  return(tbl)
  
}



##-------------------------------------------------------## Run Chromosome Compaction Analysis Function


rat.1mbp.inf  <- mwa.cont.ratio(distal.1 = 1000000, distal.2 = 1000000000, local.1 = 50000, local.2 = 1000000,
                                inc = 1000000, window = 5000000, dir = dir.contacts, exps = c("WT_HCT116_CTRL", "Pol1_6hrs_Aux", "Pol2_6hrs_Aux","1hr_ActD"))

rat.1mbp.inf  <- mwa.cont.ratio(distal.1 = 1000000, distal.2 = 1000000000, local.1 = 50000, local.2 = 1000000,
                                inc = 500000, window = 1000000, dir = dir.contacts, exps = c("WT_HCT116_CTRL", "Pol1_6hrs_Aux", "Pol2_6hrs_Aux","1hr_ActD"))

rat.250kbp.1mbp  <- mwa.cont.ratio(distal.1 = 250000, distal.2 = 1000000, local.1 = 50000, local.2 = 250000,
                                   inc = 100000, window = 100000, dir = dir.contacts,exps = c("WT_HCT116_CTRL", "Pol1_6hrs_Aux", "Pol2_6hrs_Aux","1hr_ActD"))

##-------------------------------------------------------## Write out results

write.csv(rat.1mbp.inf, file.path("/home/lmc0633/Transcription Paper/", "rat.1mbp.inf.allconds.csv"), row.names=FALSE)