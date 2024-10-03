library(strawr)

##------------------------------ available STRAWR functions ------------------------------##

## Function for reading basepair resolutions from .hic file
#~ readHicBpResolutions(system.file("extdata", "test.hic", package = "strawr"))

## Function for reading chromosomes from .hic file
#~ readHicChroms(system.file("extdata", "test.hic", package = "strawr"))

## Function for reading available normalizations from .hic file
#~ readHicNormTypes(system.file("extdata", "test.hic", package = "strawr"))

##-------------------------------- Strawr Dump Loop --------------------------------------##

# Function for dumping intrachromosomal contacts
# StrawR required for this function
# 
## fast C++ implementation of dump. Not as fully featured as the Java version. Reads the .hic file,
## finds the appropriate matrix and slice of data, and outputs as data.frame in sparse upper triangular
## format. Currently only supporting matrices.

#~ norm = "NONE"  # Must be one of NONE/VC/VC_SQRT/KR.
#~ resolution = 100000 # The bin size. By default, for BP, this is one of <2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000>
#~ unit =  "BP" # BP (BasePair) or FRAG (FRAGment)
#~ matrix = "oe" # Type of matrix to output. Must be one of observed/oe/expected. observed is observed counts, oe is observed/expected counts, expected is expected counts.

#~ hic.data.frame <- strawr::straw(norm,"/path/to/file.hic", "1", "1", unit, resolution, matrix = matrix)
# 
# Directory string should look like this: /projects/b1042/BackmanLab/juicerDir/experimentName/contact_data/cond/rep
# If no cond or reps, set to NULL
## Dump interchromosal contacts

dumpIntracontacts <- function(root.path, juicerDir, exp, conds, reps, res, norm, hic_file, matrix){
  
  #chroms <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')
  chroms = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X')
  
  if (require(strawr)){
    
    cat("strawr is loaded correctly")
    
    cat(paste0("Root Directory: ", root.path, "\n"))
    cat(paste0("Juicer Directory: ", juicerDir, "\n"))
    cat(paste0("Experiment: ", exp, "\n"))
    
    if (!is.null(conds) & !is.null(reps)){
      
      cat(paste0("condition: ", conds, "\n","\n"))
      
      for (cond in 1:length(conds)){
        
        cond <- conds[cond] ## Condition labels
        
        ## Input and output paths on Quest
        input.path = paste0(root.path, juicerDir, exp,'/',cond,'/')
        cat(paste0("\n","\n","input path: ", input.path, "\n"))
        
        output.path = paste0(root.path,juicerDir,exp,'/',cond,'/intracontacts/')
        cat(paste0("output path: ",output.path, "\n","\n"))
        
        ## Create output directory if it does not exist 
        dir.create(file.path(output.path))
        
        for (rep in 1:length(reps)){
          
          for (chrom in chroms) {
            
            print(paste0("processing intrachromosomal contacts: ", chrom, " to ", chrom, " ", matrix, "contacts at ", res, " resolution"))
            
            hic.data.frame <- strawr::straw(norm, paste0(input.path, rep, hic_file), chrom, chrom, "BP", res , matrix=matrix)
            
            output.fname = paste("rep",rep,'-chr',chrom,'_chr',chrom,'-',matrix,"_",(res/1000),'Kb_',norm,'norm.csv',sep="")
            print(output.fname)
            
            #write.csv(hic.data.frame, paste0(output.path, output.fname), row.names=FALSE)
            cat(paste0("\n","\n", "Contacts have been dumped in ", output.path, " for ", exp,"\n","\n" ))
          }
        }
      }
    } else if (!is.null(reps) & is.null(conds)) {
      
      cat(paste0("reps: ", reps, "\n"))
      
      for (rep in 1:length(reps)){
        
        rep <- reps[rep] ## Condition labels
        
        ## Input and output paths on Quest
        input.path = paste0(root.path, juicerDir, exp,'/',rep,'/')
        cat(paste0("\n","\n","input path: ",input.path, "\n"))
        
        output.path = paste0(root.path,juicerDir,exp,'/',rep,'/intracontacts/')
        cat(paste0("output path: ",output.path, "\n","\n"))
        
        ## Create output directory if it does not exist 
        dir.create(file.path(output.path))
        
        for (chrom in chroms) {
          
          print(paste0("processing intrachromosomal contacts: ", chrom, " to ", chrom, " ", matrix, "contacts at ", res, " resolution"))
          
          hic.data.frame <- strawr::straw(norm, paste0(input.path, rep, hic_file), chrom, chrom, "BP", res , matrix=matrix)
          
          output.fname = paste(rep,'-chr',chrom,'_chr',chrom,'-',matrix,"_",(res/1000),'Kb_',norm,'norm.csv',sep="")
          print(output.fname)
          
          write.csv(hic.data.frame, paste0(output.path, output.fname), row.names=FALSE)
          cat(paste0("\n","\n", "Contacts have been dumped in ", output.path, " for ", exp,"\n","\n" ))
          
        }
      }
      
    } else {
      
      cat("no reps and no conds\n")
      
      input.path = paste0(root.path, juicerDir, exp,'/')
      cat(paste0("\n","\n","input path: ",input.path, "\n"))
      
      output.path = paste0(root.path,juicerDir,exp,'/intracontacts/')
      cat(paste0("output path: ", output.path, "\n","\n"))
      
      ## Create output directory if it does not exist 
      dir.create(file.path(output.path))
      
      for (chrom in chroms) {
        
        print(paste0("processing intrachromosomal contacts: ", chrom, " to ", chrom, " ", matrix, "contacts at ", res, " resolution"))
        
        hic.data.frame <- strawr::straw(norm, paste0(input.path, hic_file), chrom, chrom, "BP", res , matrix=matrix)
        
        output.fname = paste('chr',chrom,'_chr',chrom,'-',matrix,"_",(res/1000),'Kb_',norm,'norm.csv',sep="")
        print(output.fname)
        
        write.csv(hic.data.frame, paste0(output.path, output.fname), row.names=FALSE)
        cat(paste0("\n","\n", "Contacts have been dumped in ", output.path, " for ", exp,"\n","\n" ))
        
      }
    }
  } else {
    
    stop("strawr not installed. Please install strawr before continuing")
    
  }
}

dumpIntracontacts(root.path = "/projects/b1042/BackmanLab/", juicerDir="HiC2/opt/juicer/work/",
                  exp = c("112123_HiC"), conds = c("1hr_ActD",  "Pol2_6hrs_Aux","Pol1_6hrs_Aux", "WT_HCT116_CTRL"), 
                  reps = c("mega"), res=5000, norm="KR", hic_file = "inter_30.hic", matrix="oe")

