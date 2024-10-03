

## HiC Analysis Guide

First, we run the python script in the Jupiter Notebook, <kbd>Cluster Hi-C Script Generator.ipynb</kbd>. This notebook generates shell scripts for each step of High Throughput Chromatin Conformation Capture that happens on Northwestern's Quest HPC Cluster. These scripts are placed in the local directory that you have specified and can then be moved to the cluster directory where you'll perform the analysis. After the initial function in the script, the remaining functions in <kbd>Cluster Hi-C Script Generator.ipynb</kbd> will be included at the end of this markdown

### Step 1a: Build directory structure first:

```shell
cd /projects/b1042/BackmanLab/HiC2/

mkdir /projects/b1042/BackmanLab/HiC2/opt/
mkdir /projects/b1042/BackmanLab/HiC2/opt/juicer/
mkdir /projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/
mkdir /projects/b1042/BackmanLab/HiC2/opt/juicer/restriction_sites/
mkdir /projects/b1042/BackmanLab/HiC2/opt/juicer/references/
mkdir /projects/b1042/BackmanLab/HiC2/opt/juicer/work/ ## this is the working directory
mkdir /projects/b1042/BackmanLab/HiC2/opt/juicer/chrom_sizes/

```

### Step 1b: Create Directory Tree
First, it creates a directory structure using the following function and variables:

```python
##----------------------------Quest Allocation----------------------------##

allocationID = 'b1042'
partitionName = 'genomics'
email = 'Your_Email@u.northwestern.edu'

##----------------------------Directory Variables----------------------------##

# Local machein dir here the scripts and directory structer is initially generated
output_path = '/Users/Your_Name/Documents/Your_Working_Directory/'
# Location of juicer directory on cluster
juicer_dir = '~/HiC/opt/juicer'
# Location of non-Juicer post-Hi-C analysis python/R code on cluster
analysis_code_dir = '~/HiC/code_files/'
# Location of working juicer directory on cluster
juicer_work_dir = '/projects/b1042/BackmanLab/Your_Name/HiC/opt/juicer/work/'

##----------------------------Experiment Variables----------------------------##

# Directory name of the experiment you're analyzing
experiment_name = 'Name_Of_Experiment'
# A list of condition directories for each condition to be analyzed
conds = ['Cond_1','Cond_2']
# A list of the number of replicates in each condition
reps_per_cond = [1,1]
# Name of genome to align to
genome = "hg38"
# Name of genome fa (if needed)
fasta_name = "hg38.fa"
# Chromosomes to analyze contacts (depends on genome)
chroms = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X'}
# Name of txt file with restriction enzyme cuts of genome
restriction_sites = "My_Restriction_Site_File.txt"
# Name of restriction enzyme
restriction_enzyme = "DpnII"
# Either "Rep" (analyze by replicate) or "mega" (pool all replicates together, analyze by condition)
keyword = "Rep"

##----------------------------Directory Function----------------------------##

def create_dir_struct(output_path,experiment_name,conds,reps_per_cond):

    experiment_dir = output_path+experiment_name
    if not os.path.isdir(experiment_dir):
        os.mkdir(experiment_dir)

    subdirnames = ['juicer_analysis','contact_data','contact_domains','compartment_analysis']

    for subdirname in subdirnames:
        subdir = experiment_dir+'/'+subdirname+'/'

        if not os.path.isdir(subdir):
            os.mkdir(subdir)

```


### Step 2a: Prepare indices and cut site file

##### **Get genome reference and build indices:**

if using the rDNA reference: [Construction and validation of customized genomes for human and mouse ribosomal DNA mapping](https://www.jbc.org/article/S0021-9258(23)01794-5/fulltext)
Get reference genome that contains rDNA repeats from [here](https://github.com/vikramparalkar/rDNA-Mapping-Genomes). These genomes are for mappping only. For Cut site generation, use full hg38 from UCSC. 

```shell
wget https://github.com/vikramparalkar/rDNA-Mapping-Genomes/raw/main/Human_hg38-rDNA_genome_v1.0.tar.gz

```
For general mapping, use Ensembl:

```shell
VERSION=108
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```

OR UCSC: [UCSC Genome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/) and [UCSC GTF](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/) (preferred). There are 4 different GTF references to use.

```shell
wget -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz

```
***To build the BWA indices needed for alignment:***

transfer genome.fa to your reference directory

note: Must use version bwa/0.7.17 for compatibility with Juicer. Juicer.sh CPU version with throw an error about the -5 flag (-SP5M) if you use a different version. We also need the .SA file with the other indices/references. If BWA doesn't generate this using the following script, run a second script (below first script) to generate just that file. Make sure you give this plenty of memory in the allocation

Example SLURM header below:
```shell
#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --mem=60GB
#SBATCH --ntasks-per-node=16
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 24:00:00
#SBATCH --job-name=BWA_index

```

```shell
##--SLURM header here--## 

## Navigate to the Juicer/references directory
cd /projects/b1042/BackmanLab/HiC2/opt/juicer/references

## Load the default version of BWA
module load bwa/0.7.17

## Generate indices for Juicer
bwa index genome.fa

```
if the .SA files aren't generated by the first script

```shell
##--SLURM header here--## 

## Navigate to the Juicer/references directory
cd /projects/b1042/BackmanLab/HiC2/opt/juicer/references

## Load the default version of BWA
module load bwa/0.7.17

## Generate SA file
bwa bwt2sa genome.fa.bwt genome.fa.sa

```

after indexing, you should have the following files:

```shell
# after running BWA indexing
/opt/juicer/references/genome.fasta.sa
/opt/juicer/references/genome.fasta.ann
/opt/juicer/references/genome.fasta.amb
/opt/juicer/references/genome.fasta.pac
/opt/juicer/references/genome.fasta.bwt
```

rename indices to hg38 as follows: <kbd> hg38.fa  hg38.fa.amb  hg38.fa.ann  hg38.fa.bwt  hg38.fa.fai  hg38.fa.pac  hg38.fa.sa </kbd>  for juicer.sh

##### **Get site positions:**

The next thing we need for Juicer is a text file with all of restriction fragment cut start sites. For this paper, they used HindIII, a 6 BP restriction endonuclease that cuts at A^AGCTT. We want a file that carries the coordinates of everywhere that HindIII cuts in the human genome, so that we can piece the HiC reads back together in Juicer. There are two ways that I know of to generate this file.

The data format for the file looks like:

<kbd> chr1 11160 12411 12461 12686 12829 13315 13420 13566 ..</kbd> without any breaks or characters. examples can be found [here](https://bcm.app.box.com/v/juicerawsmirror/file/94739463895) at the Juicer AWS Mirror.

The first uses an R package called HiTC and some data wrangling:

``` r
library(HiTC)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)

dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/My_Projects/Pol2_HiC/Annotation_Files/RestrictionSites"
setwd(dir)


## Vignettes and manuals are here: http://bioconductor.org/packages/release/bioc/html/HiTC.html

chroms <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
            'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18',
            'chr19','chr20','chr21','chr22','chrX')

#########################################################################################
# Example from HiTC doc for generation of restriction fragment txt file
#########################################################################################

# Juicer uses the following format for txt file with restriction fragments: chr1 3002504 3005876 3006265 3014101 3015361 3017473 ..
# chr1 is chromosome, 3002504 is start of fragment, 3005876 is start of next fragment
# Juicer only uses the start fragment. It does not need the end or length of the fragment. Compared HiTC generated..
# HindIII hg19 cut site file to hg19_HindIII_new.txt on Juicer mirror (https://bcm.app.box.com/v/juicerawsmirror/file/94738967920)

## Extract chromosome levels (exclude mito chrom)
human_chr <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)[1:24]
human_chr

## GRanges of restriction fragments after HindIII digestion
resFrag <- getRestrictionFragmentsPerChromosome(resSite="AAGCTT", chromosomes=human_chr, overhangs5=1, genomePack="BSgenome.Hsapiens.UCSC.hg38")
resFrag

## Concatenates all chromosome fragments together as a single Granges object
allRF <- do.call("c",resFrag)
allRF@ranges@start ## Can get fragment cut sites here


## Give names to all the cut sites
names(allRF) <- unlist(sapply(resFrag, function(x){paste0("HIC_", seqlevels(x), "_", 1:length(x))}))
names(allRF)

## exports a bed file with restriction fragments
export(allRF, format="bed", con="HindIII_resfrag_hg38.bed")

##------------------------------Juicer Restriction Site File Generation ------------------------------##

## Read in BED file
data <- data.frame(read.table(file.path(dir, "HindIII_resfrag_hg38.bed"), sep= "\t", header= FALSE))

## Read in a single chromosome
chr1 <- data[data$V1 == 'chr1',]

## append the chrom to the list of start sites
start <- append("chr1", chr1$V2)

## Collapse the entire vector into a single string separated by a single space
start <-paste(start,collapse=" ")

## Write out the string without quotes as a text file
write.table(start, "Test.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = F, quote=F)

##--------------------------------------For All Chromosomes--------------------------------------##
chroms <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
            'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18',
            'chr19','chr20','chr21','chr22','chrX')


## Loop through all chromosomes and extract relevant sites
file <- character()
for (chr in 1:length(chroms)){

  in.chr <- chroms[chr]
  out.chr <- data[data$V1 ==in.chr,]

  sites <- append(chroms[chr], out.chr$V2)

  sites <-paste(sites,collapse=" ")

  file <- append(file, sites)

}

## Write out the string without quotes as a text file
write.table(file, "HindIII_hg38.txt", append = TRUE, sep = "\t", dec = ".",
            row.names = F, col.names = F, quote=F)


```

The second uses Juicer's own Python script, [generate_site_positions.py](https://github.com/aidenlab/juicer/blob/main/misc/generate_site_positions.py). This script takes a human genome assembly and finds all of the locations where the enzyme you select cuts.

This script should be run on Quest - Usage is <kbd:> python generate_site_positions.py <restriction enzyme> <genome> [location]</kbd>

In Quest, copy and paste the generate_site_positions.py script in the same directory where the ref genome is stored:
```shell
## Change directory to where ref genome is stored
cd /projects/b1042/BackmanLab/ref_hg38ensembl/hg38.fa

vim generate_site_positions.py
```

```shell

module load python

python generate_site_positions.py HindIII hg38 /projects/b1042/BackmanLab/ref_hg38ensembl/hg38.fa

```
The other important annotation file you need to run Juicer is a file called <kbd> hg38.chrom.sizes</kbd> that contains the size of each of the chromosomes. This file has the following format. This can be retrieved at UCSC's online repository with the genome reference:

```
chr1	248956422
chr10	133797422
chr11	135086622
chr11_KI270721v1_random	100316
chr12	133275309
chr13	114364328
chr14	107043718 ...
```

### Part 2b: Downloading and installing Juicer

Following that, you'll want to download the CPU collection of Juicer scripts. There are some issues running the SLURM collection of scripts on Quest (I later used the SLURM version, but had to make extensive edits to the script so that it would work efficiently with NU's HPC). As of 110122, the readme in the the latest release (Juicer 2.0) says to use the Juicer 1.6 release found [here](https://github.com/aidenlab/juicer/releases/tag/1.6) since 2.0 is still under development. Do not use v.2.0

```shell
cd /projects/b1042/BackmanLab/HiC2/opt/ ## Top directory

## Clone Juicer 1.6 scripts to directory
git clone https://github.com/theaidenlab/juicer.git --branch 1.6

## copy CPU scripts to scripts directory and clean up
cd /projects/b1042/BackmanLab/HiC2/opt/juicer/CPU
cp /projects/b1042/BackmanLab/HiC2/opt/juicer/CPU /projects/b1042/BackmanLab/HiC2/opt/juicer/scripts
rm -r AWS LICENSE LSF PBS  SLURM  UGER

## Add juicer tools to scripts
cd /projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/common
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar

cd /projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/
vim juicer.sh

## On line 170 in the juicer.sh script, make sure references says /hg38/hg38.fa
hg38) refSeq="${juiceDir}/references/hg38/hg38.fa";; ## This is what is originally in the script

## Don't forget to make all your scripts executable (chmod u+r+x for all permission changes if needed)
cd /projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/common
chmod +x ./*
cd /projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/
chmod +x juicer.sh

```
That's it! Juicer is installed.

### Step 3: Generate .HiC Files
Second, we generate shell scripts for each replicate that we want to analyze using <kbd>Cluster Hi-C Script Generator.ipynb</kbd>. These scripts load _BWA_ for alignment, set the working directory to wherever the FASTQ files are stored for that replicate, call flag <kbd>-D</kbd> _(directory where Juicer software is installed)_, flag <kbd>-g</kbd> _(genome)_, flag <kbd>-y</kbd> (restriction site file with all predicted cuts based on your genome), and call flag <kbd>-p</kbd> ( a file containing the sizes of each chromosome for your genome). More on usage can be found [here](https://github.com/aidenlab/juicer/wiki/Usage). <kbd>Cleanup.sh</kbd> cleans up big repetitive files and zip fastqs. Run after you are sure the pipeline ran successfully. Output is listed [here](https://github.com/aidenlab/juicer/wiki/Running-Juicer-on-a-cluster) but the most important output are the nter.hic / inter_30.hic files( The .hic files for Hi-C contacts at MAPQ > 0 and at MAPQ >= 30, respectively).

```shell

##--SLURM header here--## 

# Load necessary modules
module load bwa

# Set your working directory
cd /projects/b1042/BackmanLab/HiC/opt/juicer/work/experiment/samplename/rep/fastq/
# Perform Juicer Analysis on fastq files:
/projects/b1042/BackmanLab/HiC/opt/juicer/scripts/juicer.sh -D /projects/b1042/BackmanLab/HiC/opt/juicer -g hg38 -y /projects/b1042/BackmanLab/HiC/opt/juicer/restriction_sites/hg38_HindIII.txt -p /projects/b1042/BackmanLab/HiC/opt/juicer/chrom_sizes/hg38.chrom.sizes

# Perform cleanup on analysis files:
/projects/b1042/BackmanLab/HiC/opt/juicer/scripts/common/cleanup.sh

# Unzip merged_nodups file to perform mega analysis:
gunzip aligned/merged_nodups.txt.gz
```

### Step 4: Merge Replicates
If not analyzing replicates as individuals, we can create statistics and a hic file from a series of replicates using the <kbd>mega.sh</kbd> script. Usage is almost identical to <kbd>Juicer.sh</kbd>. More usage can be found in the [source code](https://github.com/aidenlab/juicer/blob/main/CPU/mega.sh). Set the working directory directly above your replicates (where juicer.sh was run) and make sure that there is a 'mega' directory here for <kbd>mega.sh</kbd> to deposit merged replicates.

```shell
##--SLURM header here--## 

# Set your working directory
/projects/b1042/BackmanLab/juicer/work/112123_HiC/juicer_analysis/WT_HCT116_CTRL/

# Merge hi-c analysis by replicate into mega analysis by condition:
/projects/b1042/BackmanLab/HiC/opt/juicer/scripts/common/mega.sh -g hg38 -b HindIII
gzip /projects/b1042/BackmanLab/HiC/opt/juicer/workVas_HiC_042922/juicer_analysis/OF/mega/aligned/merged

```

If using the modified SLURM version of juicer, the following SLURM submission script is an example:

```shell
##--SLURM header here--## 

# Set your working directory
cd /projects/b1042/BackmanLab/juicer/work/112123_HiC/juicer_analysis/WT_HCT116_CTRL/

# Merge hi-c analysis by replicate into mega analysis by condition:
/projects/b1042/BackmanLab/juicer/scripts/mega.sh -D /projects/b1042/BackmanLab/juicer -q genomics -l genomics-himem -A b1042 -Q 1440 -L 8640 -T 48 -g hg38 -y /projects/b1042/BackmanLab/juicer/restriction_sites/hg38_DpnII.txt -p /projects/b1042/BackmanLab/juicer/chrom_sizes/hg38.chrom.sizes

```

**Note:**
There are several additional analysis steps to gather contact probability statistics and do contact probability analysis, if wanted. They're not included here because those steps are written in R in a separate script. For those scripts, please see <kbd>Readme.md</kbd>.

### Step 5: Compartments, TADS, and loops

At this point, we have .hic contact files for our individual replicates and merged .hic contact files that we can use for downstream analyses. From here, we can generate eigens for compartments, call loops with HICCUPs, and find TADs with arrowhead (or an equivalent TAD caller).

### Compartments
Below is a sample script for submitting scripts for all conditions and replicates to find compartments, using Juicer's [eigen](https://github.com/aidenlab/juicer/wiki/Eigenvector) and [pearson](https://github.com/aidenlab/juicer/wiki/Pearsons) tools

```shell
##--SLURM header here--## 

## Directories
cond=(1hr_ActD Pol2_6hrs_Aux)
rep=(mega Rep1 Rep2)

for i in "${cond[@]}"
do

for j in "${rep[@]}"
do


juiceDir=/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC
juiceDir2=/projects/b1042/BackmanLab/juicer/work/112123_HiC
hicDir=/juicer_analysis
outDir=/compartment_analysis
hic_file_path=$juiceDir2$hicDir/${i}/${j}/aligned/inter_30.hic
out_file_path=$juiceDir$outDir/${i}/${j}

## Make directory if DNE
if [ ! -d "$out_file_path" ] ; then
    mkdir "$out_file_path"
    chmod 777 "$out_file_path"
        echo "making compartment results directory"
fi

## Variable
norm=KR
res=500000 #500kbp to 1mbp
groupname="a$(date +%s)" # unique name for jobs in this run

chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

for k in ${chroms[*]}

do
jid=`sbatch <<- EIGEN | egrep -o -e "\b[0-9]+$"
#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=12
#SBATCH -o "${out_file_path}/${i}_eigen-${jid}.out"
#SBATCH -e "${out_file_path}/${i}_eigen-${jid}.err"
#SBATCH -J "eigen_${groupname}_${jid}"

echo "analyzing chrom ${chroms[k]} for compartments on $(date)"
echo "HiC File Path: ${hic_file_path}"

# call eigen:
/projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/common/juicer_tools eigenvector ${norm} ${hic_file_path} ${chroms[k]} BP ${res} ${out_file_path}/eigen_${chroms[k]}_${res}bp_${norm}.txt

# Call Pearsons:
/projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/common/juicer_tools pearsons ${norm} ${hic_file_path} ${chroms[k]} BP ${res} ${out_file_path}/pearsons_${chroms[k]}_${res}bp_${norm}.txt
EIGEN`

echo -e "submiting job ${jid}\n"

done

done

done


echo -e "Job list: \n"
squeue -u lmc0633 -o "%A %T %j %E %R" | column -t
echo -e "All chromosomes have been analyzed"
```

### HiC Resolution
Below is a sample script to calculate resolution for all HiC files

```shell
##--SLURM header here--## 

## Directories
cond=(1hr_ActD Pol2_6hrs_Aux)
rep=(mega Rep1 Rep2)

for i in "${cond[@]}"
do

    for j in "${rep[@]}"
    do

        ## Variables
        juiceDir2=/projects/b1042/BackmanLab/juicer/work/112123_HiC
        hicDir=/juicer_analysis

        # Calculate map resolution:
        /projects/b1042/BackmanLab/HiC2/opt/juicer/misc/calculate_map_resolution.sh $juiceDir2$hicDir/${i}/${j}/aligned/merged_nodups.txt $juiceDir2$hicDir/${i}/${i}_${j}_resolution.txt

    done

done

```
### HICCUPS
Below is a sample script to find loops using Juicer's [HICCUPS](https://github.com/aidenlab/juicer/wiki/HiCCUPS). You must use the genomics-GPU node for HICCUPS.

```shell

#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --mem=100G
#SBATCH --gres=gpu:a100:4 ## This is the GPU quest has 
#SBATCH -o hiccups-%j.out
#SBATCH -e hiccups-%j.err
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -J hiccups_juicertools

## Genomics-GPU specs: 8 x 40GB Tesla A100 GPUs available on 2 nodes (four GPUs, 52 CPU cores, and 192 GB RAM on each node)
module load cuda/cuda-10.1.2-openmpi-3.1.3
module load java/jdk1.8.0_25

## Directories
cond=(1hr_ActD Pol2_6hrs_Aux)
rep=(mega Rep1 Rep2)

for i in "${cond[@]}"


    for j in "${rep[@]}"
    do

        ## Variables

        juiceDir=/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC
        juiceDir2=/projects/b1042/BackmanLab/juicer/work/112123_HiC
        hicDir=/juicer_analysis
        outDir=/loop_analysis
        genomeID=hg38

        echo "processing cond: ${i} & rep: ${j}"
        echo "inter.HiC file is at: $juiceDir2$hicDir${i}/${j}/aligned/inter.hic"
        echo "Out directory is: $juiceDir$outDir${i}/${j}/hiccups_results"

        # hiccups annotation of contact domains:
        /projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/common/juicer_hiccups.sh -j /projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/common/juicer_tools -k KR -r 5000,10000,25000 -g $genomeID -i $juiceDir2$hicDir${i}/${j}/aligned/inter_30.hic -o $juiceDir$outDir$cond/${j}/hiccups_results

    done

done

```

### Arrowhead

Arrowhead is the domain identification algorithm used in [Rao et. al., 2014](https://www.cell.com/fulltext/S0092-8674(14)01497-4). It's the most commonly used TAD caller currently. Arrowhead usage is as follows: <kbd>arrowhead (HiC file) (output_file)</kbd>. Documentation can be found [here](https://github.com/aidenlab/juicer/wiki/Arrowhead) <kbd>-r</kbd> resolution  <kbd>-k</kbd> normalization (NONE/VC/VC_SQRT/KR). The HiC files it takes for input are those generated by <kbd>Juicer.sh</kbd> or <kbd>Mega.sh</kbd>.

Contact domain files take the following form:
```
chrom1 x1 x2 chrom2 y1 y2 color corner_score Uvar Lvar Usign Lsign
```
Explanations of each field are [here](https://github.com/aidenlab/juicer/wiki/Arrowhead)

Below is a sample script to find TADs with Arrowhead.

```shell

##--SLURM header here--## 

## Directories
cond=(1hr_ActD Pol2_6hrs_Aux)
rep=(mega Rep1 Rep2)

for i in "${cond[@]}"
do

    for j in "${rep[@]}"
    do
        ## Dirs
        juiceDir=/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC
        hicDir=/projects/b1042/BackmanLab/juicer/work/112123_HiC/juicer_analysis/
        outDir=/contact_domains/arrowhead_domains
        hic_file_path=$hicDir${i}/${j}/aligned/inter_30.hic
        out_file_path=$juiceDir$outDir${i}/${j}/inter_30_contact_domains

        if [ ! -d "$out_file_path" ] ; then
            mkdir "$out_file_path"
            chmod 777 "$out_file_path"
                echo "making domain results directory"
        fi

        echo "processing cond; ${i} & rep: ${j}"

        ## Variable
        threads=24
        norm=KR
        res=5000

        # Arrowhead annotation of contact domains:
        /projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/common/juicer_tools arrowhead -k ${norm} -r ${res} --threads ${threads} ${hic_file_path} ${out_file_path}

    done

done
```


### APA
Below is a sample script to find loops using Juicer's [Aggregate Peak Analysis](https://github.com/aidenlab/juicer/wiki/APA). You must use the genomics-GPU node for APA.

```shell

#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --mem=100G
#SBATCH --gres=gpu:a100:4 ## This is the GPU quest has 
#SBATCH -o apa-%j.out
#SBATCH -e apa-%j.err
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -J APA_juicertools

module load cuda/cuda_8.0.61
module load java/jdk1.8.0_25

## Directories
cond=/Pol1_6hrs_Aux
rep=1
juiceDir=/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC
hicDir=/juicer_analysis
outDir=/loop_analysis
hic_file_path=$juiceDir$hicDir$cond/Rep${rep}/aligned/inter_30.hic
juicer_tools_path=/projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/common/juicer_tools
out_file_path=$juiceDir$outDir$cond/Rep${rep}/hiccups_results
apa_out_file_path=${out_file_path}"/apa_results"

## Variable
genomeID=hg38
norm=KR
res=5000

echo "HiC file is at: $juiceDir$hicDir$cond/Rep${rep}/aligned/inter_30.hic"
echo "Out directory is: $juiceDir$outDir$cond/Rep${rep}/hiccups_results"

${juicer_tools_path} apa -k ${norm} -r ${res} ${hic_file_path} ${out_file_path%.*}"/merged_loops.bedpe" ${apa_out_file_path}
    

```