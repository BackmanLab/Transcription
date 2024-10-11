#########################################################################################
## Contact Probability ##
# This script is for parsing contacts and calculating contact probability
# Replace dirs and conditions with your own
#########################################################################################

rm(list = ls())
gc()

# import packages
require("ggplot2")
require("dplyr")
require("tidyr")
require("magrittr")
require("viridis")
require("tidyr")
require("ggpubr")
require("ggridges")
require("ggside")
require("plyr")

#########################################################################################
# parse contacts
#########################################################################################

# Set directory variables here
dir.contacts <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data"

## Condition labels
exps <- c("WT_HCT116_CTRL", "Pol1_6hrs_Aux", "Pol2_6hrs_Aux",  "1hr_ActD")

files <- list.files(paste0(dir.contacts,"/",exps[1],"/mega/intracontacts/"), pattern=".csv", all.files=T, include.dirs = FALSE)
files

## loop for import
for (exp in 1:length(exps)){
  
  exp <- exps[exp]
  
  filepath <- paste0(dir.contacts,"/",exp,"/mega/intracontacts/")
  cat(paste0("input path: ",filepath, "\n","\n"))
  
  output.path <- paste0(dir.contacts,"/",exp,'/mega/parsed_contacts/')
  cat(paste0("output path: ",output.path, "\n","\n"))
  
  ## Create output directory if it does not exist 
  dir.create(file.path(output.path))
  
  for (file in 1:length(files)){
    infile <- files[file]
    
    name =sapply(strsplit(infile, "_"), `[`, 1)
    chrom=sapply(strsplit(name, "-"), `[`, 2)
    cat(paste('name',chrom,'\n\n'))
    
    cat(paste('appending contact files',exp,":",infile,'\n\n'))
    
    infile <- read.table(file.path(filepath, infile), sep= ",", header= T)
    infile <-data.frame(infile)
    infile <- na.omit(infile)
    
    infile <- infile %>% 
      mutate(index= abs(infile$x-infile$y)) 
    
    infile <-na.omit(infile)
    
    index= sort(unique(abs(infile$x-infile$y)))
    cat(paste0("chrom: ",dim(index), "\n","\n"))
    
    dist <- aggregate(counts ~ index, data=infile, sum)
    
    write.csv(dist, paste0(output.path,chrom,"_contacts_5kb.csv"), row.names=FALSE)
  }
  
}

#########################################################################################
# calculate contact probability globally
#########################################################################################

contactProbability <- function(res, exps, files, dir.domains){
  
  require(dplyr)
  
  ## Set index for merging
  index <- seq(0,(res*49340), by = res)  # value based on length of chrom1 at res of 5000
  
  for (exp in 1:length(exps)){
    
    ## Generate empty data frames as containers
    all.my.files <- data.frame(row.names = index, index = index)
    cont.prob <- data.frame(row.names = index, index = index)
    
    exp <- exps[exp]
    
    filepath <- paste0(dir.domains,"/",exp,"/mega/parsed_contacts/")
    cat(paste0("input path: ",filepath, "\n","\n"))
    
    output.path = paste0(dir.domains,"/",exp,'/mega/contact_prob/')
    cat(paste0("output path: ",output.path, "\n","\n"))
    
    ## Create output directory if it does not exist 
    dir.create(file.path(output.path))
    
    for (file in 1:length(files)){
      
      infile <- files[file]
      name =sapply(strsplit(infile, "_"), `[`, 1) # get chromosome names
      cat(paste0("chrom: ",name, "\n","\n"))
      
      cat(paste('appending contact files',exp,":",infile,'\n\n'))
      
      ## introduce file
      infile <- read.table(file.path(filepath, infile), sep= ",", header= T)
      infile <-data.frame(infile)
      
      ## Calculate contact probability for each chromosome
      infile$possible.contacts <- seq((res*nrow(infile)),res, by = -res) # get all possible contacts
      cat(paste('possible contacts: ',nrow(infile),'\n\n'))
      
      infile <- infile %>% mutate(contact.prob = counts/possible.contacts) # get the contact probability (counts/possible contacts)
      infile <- infile[,c("index", "contact.prob")] 
      colnames(infile) <- c("index", name) # rename columns to chrom
      
      ## Merge all contacts by index (1D genomic distance)
      all.my.files  <- merge(all.my.files , infile, by= "index", all.x = T)
      all.my.files[is.na(all.my.files)] <- 0 # replace NA values with 0s for each merge
      cat(paste('merged all chromosomes into one data frame ','\n\n'))
      
    }
    
    ## Get the sum of all contact probabilities
    all.my.files  <- all.my.files[-c(1)] # remove index from merge.df
    cont.prob$row.sums <- rowSums(all.my.files) # sum contact probabilites for each chromosome
    cat(paste('rowsum check: ',cont.prob$row.sums[1:3],'\n\n'))
    
    ## Get col sums for each chromosome to calculate each chrom contact probability
    col.sums <- colSums(all.my.files) #df2/df2[,3]
    all.my.files <- rbind(all.my.files, "col.sums" = col.sums)
    all.my.files<- mapply('/', all.my.files, all.my.files[col.sums,]) # [tail(all.my.files, n =1) # alt way
    
    ## remove negative values and NAs
    cont.prob[cont.prob$contact.prob < 0,] <- 0
    
    ## Normalize contact probabilities
    sum <- sum(cont.prob$row.sums) # get sum of all contact probabilities
    cat(paste('Sum of contact probabilities: ',sum,'\n\n'))
    cont.prob <- cont.prob %>% mutate(contact.prob = (row.sums/sum), .keep = "unused") # get normalized contact probability
    
    write.csv(cont.prob, paste0(output.path,exp,"_contact_prob_5kb.csv"), row.names=FALSE)
    
  }
  
  return(cont.prob)
  
}

# Set directory variables here
dir.domains <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data"

## Get list of files
files <- list.files(paste0(dir.domains,"/",exps[1],"/mega/parsed_contacts/"), pattern=".csv", all.files=F, include.dirs = FALSE)
files

## Run all conditions
contactProbability(res=5000, exps <- c("WT_HCT116_CTRL", "Pol1_6hrs_Aux" ), files =files, dir.domains = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data")

## Run each condition individually
wt.cont.prob <- contactProbability(res=5000, exps <- c("WT_HCT116_CTRL"), files =files, dir.domains = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data/")
pol1.cont.prob <- contactProbability(res=5000, exps <- c("Pol1_6hrs_Aux"), files =files, dir.domains = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data/")
pol2.cont.prob <- contactProbability(res=5000, exps <- c("Pol2_6hrs_Aux"), files =files, dir.domains = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data/")
actd.cont.prob <- contactProbability(res=5000, exps <- c("1hr_ActD"), files =files, dir.domains = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data/")

#########################################################################################
# Contact probability plotting
#########################################################################################

library(tidyr)
library(ggplot2)

## Prepare data for plotting
contrast<- cbind(wt.cont.prob, pol1.cont.prob)
contrast<- contrast[,-3]

contrast$index[1] <- 1
colnames(contrast) <- c("dist", "ctrl", "aux")

## Prepare data for plotting
plot.df <- contrast %>% 
  gather(cond, prob, 2:3) 

plot.df <- plot.df %>% 
  mutate(dist = log(dist, 10)) %>% 
  mutate(prob = log(prob, 10)) %>% 
  filter(dist < 8)

## filename: pol1.rcp.mycode
p <- plot.df %>%
  ggplot(aes(x=dist , y=prob)) + 
  geom_line(size = 0.5, aes(color = cond)) +theme_classic()+
  scale_y_continuous(limits = c(-5,  -0.5), breaks = c( -7, -6, -5, -4, -3, -2, -1), labels = c(7, 6, 5,4,3,2, 1))+
  scale_x_continuous(limits = c(3.69,7.5), breaks = c(3.5,4, 4.5,5, 5.5, 6,6.5,7, 7.5))+
  scale_color_manual(values=c("#798234","#D46780")) +
  labs(title="Contact Probability: 6hrs Pol1 Degradation")+ theme(plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  xlab("Log 10 Base Pairs") + ylab("-Log10 Contact Probability ")
p

#########################################################################################
# Logarithmic Binning Contact probability
#########################################################################################

## Use this function for calculating contact probabilty with log smoothing
## The tails of contact data over long distances become irregular as data becomes
## sparser and need smoothing to interpret.

f <- function(x){
  
  return(mmand::gaussianSmooth(x, sigma=10000))
}

contactProbBin <- function(exps, dir, files){
  
  ## loop for calculating sliding window log2 ratio
  for (exp in 1:length(exps)){
    
    tbl<- data.frame()
    
    ## Loop through conditions
    exp <- exps[exp]
    
    files <- list.files(paste0(dir,"/",exp,"/mega/contact_prob"), pattern=".csv", all.files=F, include.dirs = FALSE)
    
    for (file in 1:length(files)){
      
      ## Loop through contact files in each condition
      infile <- files[file]
      
      filepath <-paste0(dir,exp,'/mega/contact_prob/')
      cat(paste0("input file path: ",filepath," for experiment ", exp, "\n","\n"))
      
      output.path = paste0(dir,exp,'/mega/scaling_analysis/')
      cat(paste0("output path: ",output.path, "\n","\n"))
      
      ## Create output directory if it does not exist 
      #dir.create(file.path(output.path))
      
      infile <- read.table(file.path(paste0(filepath), infile), sep= ",", header= T)
      infile$index<- ifelse(infile$index== 0 , infile$index[1] <- 5000 , infile$index)
      
      #infile = infile[log10(infile$index) > start & log10(infile$index) < end,]
      
      ## Generate list spaced by increment and length of infile left contacts
      
      list <- log10(unique(as.integer(pracma::logspace(log10(min(infile$index)), log10(max(infile$index)), n = 70))))
      print(list)
      
      start <- list[c(TRUE, FALSE)]
      end <- list[c(FALSE, TRUE)]
      
      cat(paste0(" Window from ",start, " - ",end, " BP", "\n","\n"))
      
      for (i in 1:length(start)) {
        
        print(start[i])
        print(end[i])
        
        test <- infile[log10(infile$index) >= (start[i]) & log10(infile$index) <= (end[i]),]
        if (nrow(test) < 1) {
          
          next
          
        } else {
          
          fitting.range <- infile[log10(infile$index) >= (start[i]) & log10(infile$index) <= (end[i]),]
          
        }
        
        print(fitting.range[1:5,])
        
        wind <- paste0("window: ", exp, " ",start[i], " - ",end[i], " BP") ## Window string
        cat(paste0("window: ", exp, " ",start[i], " - ",end[i], " BP", "\n","\n"))
        
        fitting.range <- fitting.range %>% 
          mutate(dist = f(log(index, 10))) %>% 
          mutate(prob = f(log(contact.prob, 10))) #%>% 
        #mutate(prob = log(contact.prob, 10))
        
        df <- data.frame(exp, start[i], end[i], wind, fitting.range$dist, fitting.range$prob)
        
        tbl <- rbind(tbl,df)
        
      }
      
      cat(paste(exp," is done ",'\n\n'))
      
    }
    
    colnames(tbl) <- c("cond", "start", "end", "window", "dist", "prob")
    #write.csv(tbl, paste0(output.path,exp,"_scaling",log10(inc),"inc-",log10(window), "window.csv"), row.names=FALSE)
    
    return(tbl)
  }
  
}

## Prepare data for plotting
wt.lgbin.cp <- contactProbBin( exps <- c("WT_HCT116_CTRL"), dir = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data/", files="WT_HCT116_CTRL_contact_prob_5kb.csv") 
pol1.lgbin.cp <- contactProbBin(exps <- c("Pol1_6hrs_Aux"), dir = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data/", files="Pol1_6hrs_Aux_contact_prob_5kb.csv") 
pol2.lgbin.cp <- contactProbBin(exps <- c("Pol2_6hrs_Aux"), dir = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data/", files="Pol2_6hrs_Aux_contact_prob_5kb.csv") 
actd.lgbin.cp <- contactProbBin(exps <- c("1hr_ActD"), dir = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/contact_data/", files="Pol2_6hrs_Aux_contact_prob_5kb.csv") 

data <- rbind(wt.lgbin.cp, actd.lgbin.cp )

data[sapply(data, is.infinite)] <- NA
data <- na.omit(data )

## plot log-binned contact probability # filename: pol1.smoothed.rcp.mycode
p <- data %>%
  ggplot(aes(x=dist, y=prob)) + 
  geom_line(size = 0.5, aes(color = cond)) +theme_classic()+
  #scale_y_continuous(limits = c(-6,  -1), breaks = c(-6, -5, -4, -3, -2, -1), labels = c(6, 5,4,3,2, 1))+
  #scale_x_continuous(limits = c(3.5,8.25), breaks = c(4, 4.5,5, 5.5, 6,6.5,7,7.5,8))+
  scale_color_manual(values=c("#798234","#D46780")) +theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  labs(title="Contact Probability: 1 Hr Act D 5ug/mL")+ xlab("Log 10 Base Pairs") + ylab("-Log 10 Contact Probability ")
p

