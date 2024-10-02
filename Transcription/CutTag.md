## Cut&Tag Data Processing Guide

### Intro:

Cleavage Under Targets & Tagmentation (CUT&Tag)is an ChIP-like profiling strategy in which primary antibodies are bound to chromatin proteins in situ after nuclei have been permeabilized. Secondary antibodies covalently bound to the cut-and-paste transposase Tn5 are then bound to the primary antibody. Activation of the transposase simultaneously digests DNA and adds high throughput sequencing adapters (referred to as ‘tagmentation’) for paired-end DNA sequencing. Cut&Tag libraries have a much higher signal-to-noise ratio than traditional ChIP libraries and need far fewer cells as well.

Below is documentation for set-up, background on this project, and a single script **CutTag.sh** for parsing Cut&Tag data based off of the data analysis protocol presented [here](https://yezhengstat.github.io/CUTTag_tutorial/#I_Introduction) from the Cut&Tag [protocol](https://www.nature.com/articles/s41596-020-0373-x) and [paper](https://www.nature.com/articles/s41467-019-09982-5). Additionally, there is an [nf-core pipeline](https://nf-co.re/cutandrun/3.2.2/) for Cut&Tag. Because I want  more control over the processing of the data, I will write my own wrapper based on the Henikoff Lab's data analysis protocol.

Here are related papers and resources for Cut&Tag:
- [A 2023 review on Cut&Tag](https://www.tandfonline.com/doi/full/10.1080/15592294.2023.2293411#abstract)
- [Integrative Analysis of CUT&Tag and RNA-Seq Data Through Bioinformatics: A Unified Workflow for Enhanced Insights](https://link.springer.com/protocol/10.1007/978-1-0716-4071-5_13)
- [CUT&RUNTools 2.0: a pipeline for single-cell and bulk-level CUT&RUN and CUT&Tag data analysis](https://academic.oup.com/bioinformatics/article/38/1/252/6318389)
- [Spatial-CUT&Tag: Spatially resolved chromatin modification profiling at the cellular level](https://www.science.org/doi/full/10.1126/science.abg7216?casa_token=xLnWrQG_IG0AAAAA%3AdHWxN4ylifpH74stpHLjaS4jEISVEni-8akaksmJ-kuijNiCDKb515VaHKQWkmHSk3O_kX4ZZTTh4I4)
- [Active Motif's Spike in Strategy](https://www.activemotif.com/documents/AACR-2024-Spike-In-Poster.pdf)
- [GoPeaks: histone modification peak calling for CUT&Tag](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02707-w)
- [CUT&Tag recovers up to half of ENCODE ChIP-seq peaks in modifications of H3K27](https://www.biorxiv.org/content/10.1101/2022.03.30.486382v2.full)
- [ChIPseqSpikeInFree: a package for Spike-in Free analysis](https://github.com/stjude/ChIPseqSpikeInFree)
- [Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4)
- [Another Henikoff tutorial on Cut&Tag Analysis](https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-5jyl8py98g2w/v2) 
- [Really good resource on Bash operators](https://tldp.org/LDP/abs/html/comparison-ops.html)


I submitted three groups of Cut&Tag biological replicates on 6.03.2024,  07.01.2024, and 7.19.2024. The later two submissions have technical replicates in them that should be combined. The sample identifiers used for this study are below:

**6.03.2024 - Submission 1:**

|Cells          | AB	 | Index     | NuSeq ID |
|---------------|--------|-----------|----------|
| WT	        | Pol1AB | i71-i51	 |BM101    |
| WT	        |K4Me3AB |	i71-i54	 |BM103     |
| WT	        |Pol2AB	 |i72-i51	 |BM105     |
| WT	        |K9Me3AB |i72-i53	 |BM107     |

**07.01.2024 - Submission 2:**

|Cells  Rep1    | AB	 | Index     | NuSeq ID |
|---------------|--------|-----------|----------|
|WT	            |Pol1AB  |i71-i51	|BM109 |
|WT	            |Pol2AB	 |i71-i54	|BM112 |
|WT	| Pol1AB	| i73-i51	| BM117 |
|WT	| Pol2AB	|	i73-i54	| BM120 | 


**07.19.2024 - Submission 3:**

|Cells  Rep1    | AB	 | Index     | NuSeq ID |
|---------------|--------|-----------|----------|
|WT	|Pol1AB	|i71-i51	| BM125 |
| WT	| K9AB	| i73-i51	| BM133 |
| WT	| K9AB	| i73-i54	| BM136 |

### Alignment:

Cut&Tag uses [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for alignment.

Build Bowtie indices for reference genome prior to alignment with Bowtie using the following command <kbd> bowtie2-build path_to_reference_genome.fa prefix_to_name_indexes </kbd> using default bowtie2 on quest, bowtie2/2.4.5

For mapping to rDNA repeats + human genome, use the [following reference](https://github.com/vikramparalkar/rDNA-Mapping-Genomes/blob/main/Human_hg38-rDNA_genome_v1.0_annotation.tar.gz). For general mapping, use Ensembl 

```shell
VERSION=108
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/homo_sapiens/Homo_sapiens.GRCh38.$VERSION.gtf.gz
```

OR UCSC: [UCSC Genome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/) and [UCSC GTF](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/) (preferred). There are 4 different GTF references to use.

```shell
wget -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
wget -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz

```

Build the Bowtie indices using the below script example. 

**example SLURM header for Quest HPC**

```shell
#!/bin/bash
#SBATCH -A p32171 ##--## SLURM start
#SBATCH -p genhimem
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 24:00:00
#SBATCH --job-name=build_genome
#SBATCH --output=outlog
#SBATCH --error=errlog
```


```shell
##-- header here --##

cd /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref

module load bowtie2/2.4.5

bowtie2-build /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref/hg38.fa hg38
```

GTF file can be found at ``` /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/gtf ``` and genome index can be found at ``` /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref ```

### Peak calling

There are a few options for calling peaks on Cut&Tag data. The more leniant and traditional option is [MACS2](https://github.com/macs3-project/MACS). Using MACS2 will increase the false positive rate but will capture more peaks overall. [SEACR](https://github.com/FredHutch/SEACR/tree/master) was designed with Cut&Run/TAG data (high SN) in mind. More recent and less used, [GoPeaks](https://github.com/maxsonBraunLab/gopeaks) is also specifically designed for Cut&Run/TAG data. There are several other peak callers, but for my purposes, I will use MACS2 and SEACR. [Here](https://help.pluto.bio/en/articles/choosing-between-macs2-and-seacr-for-peak-calling-with-dna-sequencing-datasets) is a short explanation on the difference between SEACR and MACS2. This [HBC workshop](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/tree/main) covers calling peaks on Cut&Run/TAG data. 

For SEACR, we need to retrieve two scripts and put them in a directory for use later:

```shell
mkdir SEACR ## in bin or scripts
wget https://github.com/FredHutch/SEACR/raw/master/SEACR_1.3.sh
wget https://github.com/FredHutch/SEACR/raw/master/SEACR_1.3.R

```
**SEACR usage and example:** 

```shell 
## usage
bash SEACR_1.3.sh experimental bedgraph [control bedgraph | numeric threshold] ["norm" | "non"] ["relaxed" | "stringent"] output prefix 

## example
bash SEACR_1.3.sh target.bedgraph 0.01 non stringent output
```

**For MACS2 and SEACR usage in this study, consult CutTag.sh script**

**Example script to call the CutTag.sh data processing script is below**

```shell

##-- header here --##

echo -e "Usage: sh $0 -f <Forward Read> -w <Path to Working Directory> -g <Path to Bowtie Indices> -r <Path to Annotation File> -c <chromosome size file name > -t <threads> -m <min MAPQ score> -s <Spike in score (optional)> -p < path to SEACR (optional)> -b < sets MAC2 to --broad (optional)> -d <remove duplicates with Picard (optional)> \n"

## with MAC2
echo "Processing data"
for file in *R1_001.fastq.gz;
do echo "${file}";
   bash CutTag.sh -f ${file} -w /projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S3/Backman26_7.19.2024/ -g /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref/hg38 -r /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/gtf/hg38.ncbiRefSeq.gtf -c hg38.chrom.sizes -t 16 -m 2 -d ;
done

## with SEACR
echo "Processing data"
for file in *R1_001.fastq.gz;
do echo "${file}";
   bash CutTag.sh -f ${file} -w /projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/test/ -g /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref/hg38 -r /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/gtf/hg38.ncbiRefSeq.gtf -p /home/lmc0633/executables/SEACR/SEACR_1.3.sh -c hg38.chrom.sizes -t 16 -m 2;
done

```
### Post-data processing

For the correlation plot, I adapted this code for use as follows:

```R
##== R command ==##

library(corrplot)

list.files('/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S1/BED/', pattern = ".bin.bed")

n=1 # c(1,2,3)
path <- '/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/'
subdir <- c("S1", "S2", "S3")
dir <- paste0(path, subdir[n], "/BED/")
hists <- list.files(dir, pattern = ".bin.bed")

reprod = c()
fragCount = NULL
for(hist in hists){
  
  if(is.null(fragCount)){
    
    fragCount = read.table(paste0(dir, hist), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
    
  }else{
    
    fragCountTmp = read.table(paste0(dir, hist), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    
  }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 

corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 0.5, cl.cex = 0.5, addCoef.col = "black", number.digits = 2, number.cex = 0.5, col = colorRampPalette(c("midnightblue","white","darkred"))(100))


```

### Downsampling reads 

I downsampled each sample to the common lowest sample because some are over-sequenced and some are undersequenced. I worry this will contribute to variability between conditions, since some samples did experience overamplification. I tried two approaches: downsampling FASTQs and downsampling BAMs where the duplicates had been removed. For the first task, I used [seqtk](https://github.com/lh3/seqtk)

usage: seqtk sample -s <seed> <path/in/fastq.gz> <#ofsamples> > <path/out/fastq.gz>

First thing to do is count the total number of reads per file:

```shell
##-- header here --##

for file in *R1_001.fastq.gz;
do echo "${file}";
READS=$(expr $(zcat ${file} | wc -l) / 4)
echo -e "\t file: $file has $READS reads  "
echo -e "\t file: $file reads: $READS " >> reads.txt
done
```

Total number of reads per sample below:

S1:

|File                                |      Reads    	   | Name      |
|------------------------------------|-------------------|-----------|
| file: K4AB-WT_S3_R1_001.fastq.gz   | reads: 68,741,393 | K4AB_WT   |
| file: K9AB-WT_S7_R1_001.fastq.gz   | reads: 2,249,616  | K9AB_WT   |
| file: P1AB-WT_S1_R1_001.fastq.gz   | reads: 4,246,499  | P1AB_WT   |
| file: P2AB-WT_S5_R1_001.fastq.gz   | reads: 21,233,193 | P2AB_WT   |

S2:

|File                              |      Reads      | Name      |
|----------------------------------|-----------------|-----------|
|	 file: P1AB_WT_S1_R1_001.fastq.gz  |reads: 2,831,267 | P1AB_WT   |
|	 file: P2AB_WT_S4_R1_001.fastq.gz  |reads: 2,904,866 | P2AB_WT   |
|	 file: P1AB_WT_S9_R1_001.fastq.gz  |reads: 408,680   | P1AB_WT   |
|	 file: P2AB_WT_S12_R1_001.fastq.gz |reads: 1,265,939 | P2AB_WT   |

S3:

|File                                                |  Reads           | Name      |
|----------------------------------------------------|------------------|-----------|
file: P1AB_WT_S3_R2_001.fastq.gz                    | reads: 45,823,405  | P1AB_WT 
file: K9AB_WT_S11_R2_001.fastq.gz                   | reads: 63,835,066  | K9AB_WT
file: K9AB_WT_S14_R2_001.fastq.gz                   | reads: 49,808,186  | K9AB_WT

Then we can normalize FASTQs based on read # using the following script 

```shell

##-- header here --##

cutoff=1000000
norm=2000000
highnorm=10000000
seed=100

## Paths
inpath=$1 ## usage sbatch <sample.sh> </inpath/>
outpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/sampled/fastq/

module load seqtk/Dec17

for file in *R1_001.fastq.gz;
do echo "${file}";

## File names
reads=$(expr $(zcat ${file} | wc -l) / 4)
R2=$(echo $file | sed 's/R1_001.fastq.gz/R2_001.fastq.gz/g')
R1out=$(echo $file | sed 's/R1_001.fastq.gz/sampled_R1_001.fastq/g')
R2out=$(echo $R2 | sed 's/R2_001.fastq.gz/sampled_R2_001.fastq/g')

if [[ $reads -lt $cutoff ]] ## If there is no path to seacr software, use MACS2
then

    echo -e "\t Sample ${file} has too few reads: ${reads} \n"

elif [[ $reads -lt $norm ]] && [[ $reads -gt $cutoff ]]
then

echo -e "\t Sample ${file} with ${reads} reads: sampled to 1,500,000 "
echo -e "\t seqtk sample -s $seed ${inpath}${file} 1500000 > ${outpath}${R1out} \n"

seqtk sample -s $seed ${inpath}${file} 1500000 > ${outpath}${R1out}
seqtk sample -s $seed ${inpath}${R2} 1500000 > ${outpath}${R2out}

elif [[ $reads -gt $highnorm ]] 
then

echo -e "\t Sample ${file} with ${reads} reads: sampled to $highnorm "
echo -e "\t seqtk sample -s $seed ${inpath}${file} $highnorm > ${outpath}${R1out} \n"

seqtk sample -s $seed ${inpath}${file} $highnorm > ${outpath}${R1out}
seqtk sample -s $seed ${inpath}${R2} $highnorm > ${outpath}${R2out}

else

echo -e "\t Sample ${file} with ${reads} reads: sampled to $norm "
echo -e "\t seqtk sample -s $seed ${inpath}${file} $norm > ${outpath}${R1out} \n"

seqtk sample -s $seed ${inpath}${file} $norm > ${outpath}${R1out}
seqtk sample -s $seed ${inpath}${R2} $norm > ${outpath}${R2out}

fi
done

gzip *.fastq


```

I checked the duplication rate using the following code and found that the rates for each sample were quite high (consistent with overamplication by PCR (my fault) and oversequencing by the core (NuSeq's fault)). [This publication](https://www.biorxiv.org/content/10.1101/2022.03.30.486382v1.full.pdf) addresses these issues and found similar duplication issues. Here is the code I used:

```R
## Summarize the duplication information from the picard summary outputs.

n=3
dir.sam <- "/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/"
subdirs <- c("S1", "S2", "S3")
sam.path <- paste0(dir.sam, subdirs[n],"/SAM/")
hists <- list.files(sam.path, pattern = ".dupMark.txt")

dupResult = c()
for(hist in hists){
  dupRes = read.table(paste0(sam.path,hist), header = TRUE, fill = TRUE)
  
  histInfo = strsplit(hist, "_")[[1]]
  dupResult = data.frame(Histone = histInfo[1], Replicate = histInfo[2], MappedFragNum_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric, DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% as.numeric * 100, EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>% mutate(UniqueFragNum = MappedFragNum_hg38 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}

```
This method for counting the number of reads in BAMs after duplicates have been removed and indexed is very fast and can be run on the command line without a submission script

```shell
## count the number of reads in indexed BAMs

for file in *.ind.bam;
do echo "${file}";
READS=$(expr $(samtools view -c ${file} ))
echo -e "\t file: $file has $READS reads  "
echo -e "\t $file \t $READS " >> reads.txt
done

```
for downsampling and normalizing dedupped BAMs, an example script is below

```shell

##-- header here --##

## script for downsampling BAMs
cd /projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S3/BAM/
module load samtools

## paths and variables
outdir=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/sampled/BAM
seed=100

reads=3000000 ## K9 S3
array=( K9AB_WT_S11.ind.bam  K9AB_WT_S14.ind.bam ) 

## function
for bam in "${!array[@]}"; do

bamS=$(echo "${array[bam]}" | sed 's/.ind.bam/.S3.samp.ind.bam/g')
fraction=$(samtools idxstats "${array[bam]}" | cut -f3 | awk -v ct=$reads 'BEGIN {total=0} {total += $1} END {print ct/total}')
samtools view -b -s ${seed}$(echo $fraction | sed 's/^0*//') "${array[bam]}" > ${outdir}/$bamS

echo -e "\n  samtools view -b -s ${seed}$(echo $fraction | sed 's/^0*//') "${array[bam]}" > ${outdir}/$bamS"
echo -e "sampling ${array[bam]} at ${fraction} to ${bamS} \n"

nozero=$(echo $fraction | sed 's/^0*//')
done

```

### merging sampled bams to call single peak set for each

Once my BAMs were downsampled, now that BAMs have been sampled, I merged them so that I could maximize signal for peak calling. The following script is an example:

```shell

##-- header here --##

in=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/sampled/BAM/
out=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/merged/BAM/

suffix1=.samp.ind.bam
suffix2=.mrg.bam
cores=24

#Merge input BAMS
module load samtools

samtools merge ${out}P1AB_WT${suffix2} ${in}P1AB_WT_S1${suffix1} ${in}P1AB_WT_S1.S1${suffix1} ${in}P1AB_WT_S1.S2${suffix1} ${in}P1AB_WT_S3.S3${suffix1} -@ $cores

```

### Generating coverage profiles

I found the overlap between genes annotated in the Ensembl GTF in and outside of NADs using Bedtools intersect, and used these to compute the enrichment profiles for K9, K4, Pol 1, and Pol 2 bigwigs in and outside of NADs at TSS of genes.

```shell

## Paths to annotations and nads
path1=/projects/b1042/BackmanLab/Lucas/090124_CnT/genome/refR/Human_hg38-rDNA_genome_v1.0_annotation/
path2=/projects/b1042/BackmanLab/Lucas/090124_CnT/NADs/

## get GTF intersection using bed coordinates
sed 's/^chr\|%$//g' ${path2}4DNFI4HQPGVC.bed  > 4DNFI4HQPGVC.bed

bedtools intersect -wa -a ${path1}Homo_sapiens.GRCh38.111.gtf -b 4DNFI4HQPGVC.bed > gene.nad.ov.gtf
bedtools intersect -wa -a ${path1}Homo_sapiens.GRCh38.111.gtf -b 4DNFI4HQPGVC.bed -v > gene.nad.nonov.gtf
####
```

```shell

##-- header here --##

module load deeptools

##== linux command ==##
outpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/merged/heatmaps/
bgpath=/projects/p32171/101823_rnaseq/Backman19_10.18.2023/totalRNA_coverage/
ndpath=/projects/p32171/101823_rnaseq/ref_hg38ensembl/

cores=10

computeMatrix scale-regions  -S ${bgpath}wt.merge.rpgc.bigWig ${bgpath}actd.merge.rpgc.bigWig  ${bgpath}pol1.merge.rpgc.bigWig ${bgpath}pol2.merge.rpgc.bigWig \
               -R ${ndpath}gene.nad.ov.gtf ${ndpath}gene.nad.nonov.gtf  \
              --skipZeros -o ${outpath}rnaseq_nads_mat.gz -p $cores -a 2000 -b 2000 

plotHeatmap -m ${outpath}rnaseq_nads_mat.gz  -out ${outpath}rnaseq_nads_heatmap.png --colorMap viridis --sortUsing sum --legendLocation 'lower-center' --regionsLabel "NADs" "non-NADS" --startLabel "Start" --endLabel "End" --samplesLabel "WT" "ActD" "Pol1" "Pol2"

```

I found the overlap between Pol 1 MACS2 peaks in and outside of NADs using Bedtools intersect (similar to above). I used these to compute the enrichment profiles for Pol1 peaks using K9, K4, Pol 1, and Pol 2 bigwigs in and outside of NADs

```shell

##-- header here --##

module load deeptools

##== linux command ==##
inpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S2/peaks/
outpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/merged/heatmaps/
bgpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/merged/coverage/
ndpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/merged/nad.beds/

cores=10

computeMatrix reference-point -S ${bgpath}P1AB_WT.cpm.bigWig ${bgpath}P2AB_WT.cpm.bigWig ${bgpath}K4AB_WT.cpm.bigWig ${bgpath}K9AB_WT.cpm.bigWig \
               -R ${ndpath}pol1.nads.bed ${ndpath}pol1.nonnads.bed  \
              --skipZeros -o ${outpath}k4nads_mat.gz -p $cores -a 500 -b 500 --referencePoint center

plotHeatmap -m ${outpath}pol1.nads_mat.gz -out ${outpath}pol1.nads_mat_heatmap.png --colorMap viridis --sortUsing sum --legendLocation 'lower-center' --regionsLabel "NADs" "non-NADS" --startLabel "Peak Start" --endLabel "Peak End" --samplesLabel "P1AB" "P2AB" "K4AB" "K9AB"

```
