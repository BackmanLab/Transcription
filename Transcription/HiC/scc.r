#########################################################################################
# HiC rep
#########################################################################################

library(hicrep)
library(readr)
library(pheatmap)

rm(list = ls())
gc()

## directories and variables 
dir.1 <- "/projects/p32171/HiC2/opt/juicer/work/112123_HiC/juicer_analysis/"
dir.2 <- "/projects/p32171/juicer/work/112123_HiC/juicer_analysis/juicer_analysis/"
list.files(paste0(dir.1), all.files=F, include.dirs = FALSE)
hic ="inter_30.hic"

exps <- c("1hr_ActD",  "Pol2_6hrs_Aux","Pol1_6hrs_Aux", "WT_HCT116_CTRL")

## -- introduce data -- ##

## List of reps
wt.rep1 <- paste0(dir.2,exps[4],"/Rep1/aligned/",hic)
wt.rep2<- paste0(dir.2,exps[4],"/Rep2/aligned/",hic)
pol1.rep1 <- paste0(dir.1,exps[3],"/Rep1/aligned/",hic)
pol1.rep2 <- paste0(dir.1,exps[3],"/Rep2/aligned/",hic)
actd.rep1 <- paste0(dir.2,exps[1],"/Rep1/aligned/",hic)
actd.rep2 <- paste0(dir.2,exps[1],"/Rep2/aligned/",hic)
pol2.rep1 <- paste0(dir.2,exps[2],"/Rep1/aligned/",hic)
pol2.rep2 <- paste0(dir.2,exps[2],"/Rep2/aligned/",hic)

##-------------------------------------------------------## Compute all chromosomes for all conditions

# Sample file paths
array <- c(wt.rep1, wt.rep2, pol1.rep1, pol1.rep2, actd.rep1, actd.rep2, pol2.rep1, pol2.rep2)
# Chromosomes
chromosomes <- paste0(c(as.character(1:22), "X"))
# Resolution
res <- 25000
# Path with data
dir_data <- "/home/lmc0633/Transcription Paper"
# File name to save the results
fileNameOut1 <- file.path(dir_data, paste0("HiCrep_Matrix_", res, ".csv"))
# Empty matrix to store HiCrep results for each sample. Filled with 1s initially
mtx_corr <- matrix(data = 1, nrow = length(array), ncol = length(array))
# Process each pair of samples
for (i in 1:length(array)) {
	for (j in 1:length(array)) {
	  # Need to process one index pair, the reverse comparison is identical
	  if (i > j) {
	    print(paste(array[i], array[j]))
	    fileNameIn1 <- array[i]
	    fileNameIn2 <- array[j]
	    # Empty vector to store chromosome-specific stratified correlation coefficients
	    chromosomes_corr <- vector(mode = "numeric", length = length(chromosomes))
	    for (k in 1:length(chromosomes)) {
	      mat1.chr <- hic2mat(fileNameIn1, chromosome1 = chromosomes[k], chromosome2 = chromosomes[k], resol = res, method = "NONE") 
	      mat2.chr <- hic2mat(fileNameIn2, chromosome1 = chromosomes[k], chromosome2 = chromosomes[k], resol = res, method = "NONE") 
	      # Depending on matrix dimensions, compare them in different order. Workaround of https://github.com/TaoYang-dev/hicrep/issues/70
	      if (nrow(mat1.chr) > nrow(mat2.chr)) {
	        scc.out = get.scc(mat2.chr, mat1.chr, resol = res, h = 5, lbr = 0, ubr = 5000000)
	      } else {
	        scc.out = get.scc(mat1.chr, mat2.chr, resol = res, h = 5, lbr = 0, ubr = 5000000)
	      }
	      # Collect chromosome-specific SCC
	      chromosomes_corr[k] <- scc.out$scc
	    }
	    # Average chromosome-specific SCCs for a sample. Store symmetrically
	    mtx_corr[i, j] <- mtx_corr[j, i] <- mean(chromosomes_corr)
	  }
	}
}

# Add sample names as rows and columns
colnames(mtx_corr) <- rownames(mtx_corr) <- array
# Save the data
write.csv(as.data.frame(mtx_corr), fileNameOut1)

names <- c("WT Rep1", "WT Rep2", "Pol1 Rep1", "Pol1 Rep2", "ActD Rep1", "ActD Rep2", "Pol2 Rep1", "Pol2 Rep2")
mtx_corr <- data.frame(mtx_corr)
colnames(mtx_corr) <- rownames(mtx_corr) <- names

## Color scales
colors <- colorRampPalette(c("#D46780","white","#798234"))(64)

# plot heat map ## sample.corr.plot.nopol2 ## sample.corr.plot.pol1only
p<-pheatmap(mtx_corr, 
            color = colors,
            fontsize = 16,
            fontsize_row = 16, 
            fontsize_col = 16,
            cellwidth = 40,
            cellheight = 40,
            cluster_rows = T,
            cluster_cols=T,
            filename = "/home/lmc0633/Transcription Paper/scc.corr.tiff",
            main = "Stratum-adjusted Correlation Coefficient: All Samples",
            annotation_legend = TRUE)
p


##-------------------------------------------------------## Get all chroms

# Sample file paths
array <- c(wt.rep1, wt.rep2, pol1.rep1, pol1.rep2, actd.rep1, actd.rep2, pol2.rep1, pol2.rep2)
# Chromosomes
chromosomes <- paste0(c(as.character(1:22), "X"))
# Resolution
res <- 25000
# Path with data
dir_data <- "/home/lmc0633/Transcription Paper"
# File name to save the results
fileNameOut2 <- file.path(dir_data, paste0("HiCrep_AllChroms_", res, ".csv"))
# Empty DF for storage
corr.df <- data.frame()
# Process each pair of samples
for (i in 1:length(array)) {
  for (j in 1:length(array)) {
    # Need to process one index pair, the reverse comparison is identical
    if (i > j) {
      print(paste(array[i], array[j]))
      fileNameIn1 <- array[i]
      fileNameIn2 <- array[j]
      # Empty vector to store chromosome-specific stratified correlation coefficients
      chromosomes_corr <- vector(mode = "numeric", length = length(chromosomes))
      for (k in 1:length(chromosomes)) {
        mat1.chr <- hic2mat(fileNameIn1, chromosome1 = chromosomes[k], chromosome2 = chromosomes[k], resol = res, method = "NONE") 
        mat2.chr <- hic2mat(fileNameIn2, chromosome1 = chromosomes[k], chromosome2 = chromosomes[k], resol = res, method = "NONE") 
        # Depending on matrix dimensions, compare them in different order. Workaround of https://github.com/TaoYang-dev/hicrep/issues/70
        if (nrow(mat1.chr) > nrow(mat2.chr)) {
          scc.out = get.scc(mat2.chr, mat1.chr, resol = res, h = 5, lbr = 0, ubr = 5000000)
        } else {
          scc.out = get.scc(mat1.chr, mat2.chr, resol = res, h = 5, lbr = 0, ubr = 5000000)
        }
        # Collect chromosome-specific SCCs for each rep and store
        chromosomes_corr[k] <- scc.out$scc
        
        rep1 = sapply(strsplit(fileNameIn1, "/juicer_analysis/"), `[`, 2)
        rep2 = sapply(strsplit(fileNameIn2, "/juicer_analysis/"), `[`, 2)
        
        mtx.df <- data.frame(rep1 = rep1, rep2 = rep2, chr = chromosomes[k], scc = chromosomes_corr[k])
        corr.df <- rbind(corr.df, mtx.df )
      }
    }
  }
}


# Save the data
write.csv(as.data.frame(corr.df), fileNameOut2)
