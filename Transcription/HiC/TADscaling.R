## Script for quantifying contact scaling for intervals

rm(list = ls())
gc()

# import packages
library(tidyverse)
library(forcats)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# Set directory variables here
#dir.domains <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_domains/arrowhead_domains"
dir.loops <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/loop_analysis/"

## Condition labels
exps <- c("1hr_ActD", "WT_HCT116_CTRL", "Pol1_6hrs_Aux", "Pol2_6hrs_Aux")

#files <- list.files(paste0(dir.domains,"/",exps[1],"/mega/inter_30_contact_domains/"), pattern=".bedpe", all.files=T) ## TADs from our data
files <- list.files(paste0(dir.loops,"/",exps[1],"/mega/hiccups_results/"), pattern="loops.bedpe", all.files=T) # loops from our data
files

## loop for import all domains

all.my.files <- c()
for (exp in 1:length(exps)){
  
  exp <- exps[exp]
  #filepath <- paste0(dir.domains,"/",exp,"/mega/inter_30_contact_domains") # tads
  filepath <- paste0(dir.loops,"/",exp,"/mega/hiccups_results") # loops
  
  
  for (file in 1:length(files)){
    infile <- files[file]
    
    cat(paste('appending .BED files',exp,":",infile,'\n\n'))
    
    infile <- read.table(file.path(filepath, infile), sep= "\t", header= FALSE)
    
    infile <-data.frame(infile)
    #infile <- infile[,c(1:4)] # tads
    infile <- infile[,c(1, 22, 23)] # loops
    
    infile <- infile %>% 
      mutate(V1=str_trim(infile$V1, "right"))%>% # trim white space from chrom names
      mutate(name =exp) %>%
      #mutate(size= abs(infile$V2-infile$V3)) # tads
      mutate(size= abs(infile$V22-infile$V23)) # loops
    
    #infile <- infile [,-c(4)] # loops
    all.my.files <- rbind(all.my.files,infile)
    
  }
}

my.domains <- all.my.files
colnames(my.domains) <- c("chr","start","end", "name","size") # loops

# Set directory variables here
dir <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data" # our contact data

## Condition labels
exps <- c("WT_HCT116_CTRL","1hr_ActD",  "Pol1_6hrs_Aux", "Pol2_6hrs_Aux")

files <- list.files(paste0(dir,"/",exps[1],"/mega/intracontacts/"), pattern=".csv", all.files=T, include.dirs = FALSE)
files

tbl<- data.frame()
for (exp in 1:length(exps)){
  
  res = 5000
  
  exp <- exps[exp]
  
  filepath <- paste0(dir,"/",exp,"/mega/intracontacts")
  cat(paste0("input path: ",filepath, "\n","\n"))
  
  #output.path <- paste0(dir.domains,"/",exp,'/mega/parsed_contacts/')
  #cat(paste0("output path: ",output.path, "\n","\n"))
  
  ## Create output directory if it does not exist 
  #dir.create(file.path(output.path))
  
  for (file in 1:length(files)){
    infile <- files[file]
    
    name=sapply(strsplit(infile, "_"), `[`, 1)
    chrom=sapply(strsplit(name, "-"), `[`, 2)
    chrom=sapply(strsplit(chrom, "chr"), `[`, 2)
    cat(paste('chrom',chrom, '\n\n'))
    
    cat(paste('appending contact files',exp,":",infile,'\n\n'))
    
    infile <- fread(file.path(filepath, infile), sep= ",", header= T)
    infile <-data.frame(infile)
    infile <-na.omit(infile) ## Remove NAs
    print(infile)
    
    ## TAD file goes in here
    domains <- my.domains
    
    domain.test <-domains[domains$chr == chrom & domains$name == exp,] ## match contacts with correct with domains
    cat(paste('Test data index: ',domain.test[1] ,'\n\n'))
    
    if (nrow(domain.test) < 1 | all(is.na(domain.test))==T ){
      cat(paste('domains empty: ',domain.test[1],'\n\n'))
      next
    } else {
      domains <-domains[domains$chr == chrom & domains$name == exp,] ## match contacts with correct with domains
    }
    
    cat(paste('using TADs from chrom: ',domains$chr," & experiment: ",domains$name ,'\n\n'))
    
    for (i in 1:nrow(domains)) {
      
      #cat(paste('domain start: ',domains$start[i]," domain end: ",domains$end[i],'\n\n'))
      
      test <- infile[infile$x >= (domains$start[i]) & infile$x <= (domains$end[i]),]
      #cat(paste('Test data index: ',test$x[1] ,'\n\n'))
      
      if (nrow(test) < 1 | all(is.na(test))==T ){
        cat(paste('Test data skipped because test data index: ',test$x[1],'\n\n'))
        next
      } else {
        ## Collect rows between the size of increment and increment + window
        fitting.range <- infile[infile$x >= (domains$start[i]) & infile$x <= (domains$end[i]),]
        domain.size <- abs(domains$start[i]-domains$end[i]) ## Get domain size TADs
      }
      
      ## extract distances from fitting range
      fitting.range  <- fitting.range  %>% 
        mutate(index= abs(fitting.range$x-fitting.range$y)) 
      
      fitting.range <-na.omit(fitting.range)
      
      index= sort(unique(abs(fitting.range$x-fitting.range$y)))
      dist <- aggregate(counts ~ index, data=fitting.range, sum)
      
      ## Limit contacts to size of each domain
      dist <- dist[dist$index < as.numeric(domain.size),] ## TADs
      #dist <- dist[dist$index < as.numeric(3000000),] ## Loops
      
      dist$possible.contacts <- seq((res*nrow(dist)),res, by = -res) # get all possible contact
      #cat(paste('possible contacts: ',min(dist$possible.contacts)," - ", max(dist$possible.contacts),'\n\n'))
      
      dist <- dist %>% mutate(contact.prob = counts/possible.contacts)  # get the contact probability (counts/possible contacts)
      col.sum <- sum(dist$contact.prob)
      #cat(paste('colsum: ',col.sum,'\n\n'))
      
      dist <- dist %>% mutate(contact.prob = contact.prob/col.sum)
      dist$index<- ifelse(dist$index== 0 , dist$index[1] <- 1 , dist$index) ## replace 0 in index with 1 for parsing
      
      range <- dist %>% 
        mutate(distance = log(index, 10)) %>% 
        mutate(prob = log(contact.prob, 10))
      
      ## filter out any infinite or zero values
      range[is.infinite(range$prob),] <- NA
      range <- na.omit(range)
      
      ## fit simple linear regression model (OLS)
      model <- lm(prob~distance, data=range, na.action="na.exclude") ## ŷ = b0 + b1x, ŷ: The estimated response value, b0: The intercept of the regression line, b1: The slope of the regression line
      
      ## |s| is the slope, coefficient 2
      slope <- model$coefficients[2]
      intercept <- model$coefficients[1]
      
      ## Collect data 
      df <- data.frame(exp,chrom, domains$start[i],domains$end[i], log10(domains$start[i]),log10(domains$end[i]), slope, intercept, dist$contact.prob[i], dist$counts[i])
      tbl <- rbind(tbl,df)
      
      #write.csv(dist, paste0(output.path,chrom,"_contacts_5kb.csv"), row.names=FALSE)
    }
    
  }
  
}

colnames(tbl) <- c("cond", "chrom", "start", "end","log10start", "log10end", "slope", "intercept", "contact.prob", "counts")
write.csv(tbl, paste0("/home/lmc0633/","scaling.csv"), row.names=FALSE)
