#!/bin/bash
##################################################################
# Author: Lucas Carter                                                   
# Email: lucascarter2025@u.northwestern.edu                                                          
# PI: Vadim Backman                                                   
# Description: 
# This is a modified version of CutTag.sh for handling 
# downsampled BAMs. refer to CutTag.sh for full usage
################################################################

##--------------------------------------------------------------## Module load

## Customize these for your environment
module load samtools/1.6
module load java/jdk1.8.0_191
module load bedtools/2.30.0 
module load deeptools/3.1.1

##--------------------------------------------------------------## Script begin

### Usage function tells users how to run the software
helpFunction()
{
      echo "*********************************** how to use CutTag.sh ***********************************"
      echo -e "Usage: sh $0 -f <Forward Read> -w <Path to Working Directory> -g <Path to Bowtie Indices> -r <Path to Annotation File> -c <chromosome size file name > -t <threads> -m <min MAPQ score> -s <Spike in score (optional)> -p < path to SEACR (optional)> -b < sets MAC2 to --broad (optional)> -d <remove duplicates with Picard (optional)> \n"
      echo -e "Run script in directory where fastqs are located \n"
      echo "This script completes several functions: Read QC, Quality Control, Alignment, Sorting, and read count of features:"
      echo -e "\t -- Quality Control: eliminates the adaptor and low quality reads using trim_galore."
      echo -e "\t -- Mapping: Mapping via Bowtie2 for all features"
      echo -e "\t -- Duplication: Marking duplicates and removing them (optional)"
      echo -e "\t -- Sorting: Sort the bam files, mark uniquely mapping reads, generate coverage files "
      echo -e "\t -- Conversion: Convert clean BAMs to BED and BEDGRAPH format "
      echo -e "\t -- Peak Calling: Call peaks using MACS2 or SEACR"
      echo -e "\t -- Read QC: Checks quality of reads using Fastqc \n \n"
      echo -e "\t-h help \n"

      echo "For more detail information, please feel free to contact: lucascarter2025@u.northwestern.edu"
      echo "**************************************************"
   exit 1 # Exit script after printing help
}


##--------------------------------------------------------------## options

while getopts "f:w:r:c:t:s:p:b" opt
do
   case $opt in
      f) R1=$OPTARG ;; ## BAM | string
      w) wdir=$OPTARG ;; ## working directory | string
      r) gtfdir=$OPTARG ;; ## path/to/gtf_file | string
      c) csize=$OPTARG ;; ## name of chrom sizes file | string
      t) threads=$OPTARG ;; ## number of threads to HPC | integer
      s) spikein=$OPTARG ;; ## spike in value | integer
      p) seacr=$OPTARG ;; ## call peaks with SEACR instead of MACS2 | path to SEACR.sh/.r | string
      b) broad=1 ;; ## broad peaks option for MACS2 | no input
      ?) echo "Unknown argument --${OPTARG}"; helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "${R1}" ] 
then
   echo "*** error: must specify forward read ***";
   helpFunction
fi

if [ -z "${wdir}" ] 
then
   echo "*** error: working directory where FASTQs are located must be provided ***";
   helpFunction
fi

if [ -z "${gtfdir}" ] 
then
   echo "*** error: path to genome annotation file must be provided ***";
   helpFunction
fi

if [ -z "${csize}" ] 
then
   echo "*** error: must provide chromosome size file ***";
   helpFunction
fi

if [ -z "${threads}" ] 
then
   echo "*** error: number of processors must be provided ***";
   helpFunction
fi


##--------------------------------------------------------------## Paths

# Get absolute path for scripts and results; check if required scripts exist

# Set working directory by moving script to directory where sampled BAMs are located
cd ${wdir}

echo "$(date): Processing alignment of BAM files"
echo -e "Generating directories. \n"

# make sure references are prepared correctly
root="$(dirname "$(pwd)")/"

echo "Workign directory is ${root}"
echo "Genome annotation file is located at ${gtfdir} \n"

# Make coverage directory
covfold="${root}coverage"
[ ! -d $covfold ]&&mkdir $covfold

# Make .BAM directory
bamfold="${root}BAM"
[ ! -d $bamfold ]&&mkdir $bamfold

# Make .BED directory
bedfold="${root}BED"
[ ! -d $bedfold ]&&mkdir $bedfold

# Make .BEDGRAPH directory
bgfold="${root}BEDGRAPH"
[ ! -d $bgfold ]&&mkdir $bgfold

# Make peaks directory
pkfold="${root}peaks"
[ ! -d $pkfold ]&&mkdir $pkfold

# Make peaks directory
rsfold="${root}results"
[ ! -d $rsfold ]&&mkdir $rsfold

##--------------------------------------------------------------## 1. Sort and statistics

cd ${bamfold} ## move to SAM alignment folder
B2=$(echo $R1 | sed 's/.samp.ind.bam/.ind.bam/g')

###
echo "1. sorting BAM. Retaining the uniquely mapping reads for ${R1}. Generating aligmnent statistics."

samtools sort ${bamfold}/${R1} > ${bamfold}/${B2}
samtools index ${bamfold}/${B2}  ## index reads
samtools flagstat ${bamfold}/${B2} > ${bamfold}/${B2}.FLAGSTAT.txt ## get alignment scores

echo -e "samtools finished processing ${R1} to ${B2} \n"
###

##--------------------------------------------------------------## 2. unscaled BAM coverage

cd ${bamfold}
BW1=$(echo $B2 | sed 's/.ind.bam/.cpm.bigWig/g')
BW2=$(echo $B2 | sed 's/.ind.bam/.rpgc.bigWig/g')

###
echo "2. Generating coverage files for ${B2}"
echo -e "\t bamCoverage --bam ${bamfold}/${B2} --normalizeUsing CPM --outFileName ${covfold}/${BW1} --binSize 1 --numberOfProcessors max"
echo -e "\t bamCoverage --bam ${bamfold}/${B2} --normalizeUsing  RPGC --outFileName ${covfold}/${BW2} --effectiveGenomeSize 2747877702 --binSize 1 --numberOfProcessors max"

# Bam coverage to generate BigWig
bamCoverage --bam ${bamfold}/${B2} --normalizeUsing CPM --outFileName ${covfold}/${BW1} --binSize 1 --numberOfProcessors max
bamCoverage --bam ${bamfold}/${B2} --normalizeUsing  RPGC --outFileName ${covfold}/${BW2} --effectiveGenomeSize 2747877702 --binSize 1 --numberOfProcessors max

echo -e "BAM coverage files generated for ${B2} using CPM and RPGC normalization \n"
### 

##--------------------------------------------------------------## 6. Bedtools conversion

cd ${bamfold} ## move to BAM folder
BD1=$(echo $B2 | sed 's/.ind.bam/.bed/g')
BD2=$(echo $BD1 | sed 's/.bed/.filt.bed/g')
BD3=$(echo $BD2 | sed 's/.filt.bed/.fragments.bed/g')
BD4=$(echo $BD3 | sed 's/.fragments.bed/.bin.bed/g')
B2N=$(echo $B2 | sed 's/.ind.bam/.n.bam/g')

###
echo "3. Converting ${B2} to BED files ${BD1}."
echo -e "\t bedtools bamtobed -i ${bamfold}/${B2} -bedpe >${bedfold}/${BD1}"
echo -e "\t awk '$1==$4 && $6-$2 < 1000 {print $0}' ${bedfold}/${BD1} >${bedfold}/${BD2} "
echo -e "\t cut -f 1,2,6 ${bedfold}/${BD2}  | sort -k1,1 -k2,2n -k3,3n  >${bedfold}/${BD3}"
echo -e "\t awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ${bedfold}/${BD3}  | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >${bedfold}/${BD4}"
echo -e "\t samtools view -F 0x04 ${samfold}/${S1} | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > ${samfold}/${S1}_fragmentLen.txt"

samtools sort -n  ${bamfold}/${B2} > ${bamfold}/${B2N} ## Sort by read alignment flag for
bedtools bamtobed -i ${bamfold}/${B2N} -bedpe >${bedfold}/${BD1}

## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${bedfold}/${BD1} >${bedfold}/${BD2} 

## Only extract the fragment related columns
cut -f 1,2,6 ${bedfold}/${BD2}  | sort -k1,1 -k2,2n -k3,3n  >${bedfold}/${BD3}

## Bin genome for correlation calculation later
binLen=500
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ${bedfold}/${BD3}  | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >${bedfold}/${BD4}

## Extract the 9th column from the alignment sam file which is the fragment length
S1=$(echo $B2 | sed 's/.ind.bam/.sam/g')
samtools view -h ${bamfold}/${B2} -o ${bamfold}/${S1}
samtools view -F 0x04 ${bamfold}/${S1} | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > ${bamfold}/${S1}_fragmentLen.txt

echo -e "Done converting ${B2} to ${BD2} and ${BD3}. Generated binned genome ${BD4} for correlation analysis. \n"
###

##--------------------------------------------------------------## 7. Bedgraph generation

gendir=/projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref/hg38
gpath=$(echo $gendir | rev | cut -d'/' -f2- | rev) ## remove last part of directory
chromsize=${gpath}/${csize}
cd ${bedfold}

###
echo "4. Converting ${BD3} to BEDGRAPH for peak calling."
echo -e "\t chromosome size path is ${chromsize}"

if [[ -z "${spikein}" ]] ## If there is no spikein generate bedgraph without
   then

   BG1=$(echo $BD3 | sed 's/.fragments.bed/.bedgraph/g')
	echo -e "\t Generating BEDGRAPH without Spike-In"
   echo -e "\t bedtools genomecov -bg -i ${bedfold}/${BD3}  -g ${chromsize} > ${bgfold}/${BG1}"

   bedtools genomecov -bg -i ${bedfold}/${BD3}  -g ${chromsize} > ${bgfold}/${BG1}
   
   echo -e "Done converting ${BD3} to ${BG1}. \n"
   ###  

   else

   BG1=$(echo $BD3 | sed 's/.fragments.bed/.bedgraph/g')
   BW3=$(echo $BD3 | sed 's/.fragments.bed/.scaled.bigWig/g')

   libSize=`cat ${bedfold}/${BD3}|wc -l`
   scale=`echo "15000000/(${libSize}*${spikein})" | bc -l`
   echo -e "\t Generating BEDGRAPH and BIGWIG with scaled value from Spike-In: ${scale}"
   echo -e "\t bedtools genomecov -bg -scale ${scale} -i ${bedfold}/${BD3}  -g ${chromsize} > ${bgfold}/${BG1}"
   echo -e "\t bamCoverage --bam ${bamfold}/${B2} --scaleFactor ${scale} --normalizeUsing  RPGC --outFileName ${covfold}/${BW3} --effectiveGenomeSize 2747877702 --binSize 1 --numberOfProcessors max"

   bedtools genomecov -bg -scale ${scale} -i ${bedfold}/${BD3}  -g ${chromsize} > ${bgfold}/${BG1}
   bamCoverage --bam ${bamfold}/${B2} --scaleFactor ${scale} --normalizeUsing  RPGC --outFileName ${covfold}/${BW3} --effectiveGenomeSize 2747877702 --binSize 1 --numberOfProcessors max
   
   echo -e "Done converting ${BD3} to ${BG1} and ${BW3}. \n"
   ###   
   
fi

##--------------------------------------------------------------## 8. Peak calling

###

if [ -z "${seacr}" ] ## If there is no path to seacr software, use MACS2
then

cd ${bamfold}
P1=$(echo $BG1 | sed 's/.bedgraph/.macs2/g')
module load MACS2/2.2.9.1

echo "5. Calling MAC2 Peaks on ${B2}"
echo -e "\t macs2 callpeak -t ${bamfold}/${B2} -g hs -f BAMPE -n ${P1} --outdir ${pkfold} -q 0.1 --keep-dup all 2>${pkfold}/${P1}.macs2Peak_summary.txt"

macs2 callpeak -t ${bamfold}/${B2} -g hs -f BAMPE -n ${P1} --outdir ${pkfold} -q 0.1 --keep-dup all 2>${pkfold}/${P1}.macs2Peak_summary.txt

elif [[ $broad -eq 1 ]]
then

cd ${bamfold}
P1=$(echo $BG1 | sed 's/.bedgraph/.macs2/g')
module load MACS2/2.2.9.1

echo "6. Calling MAC2 Peaks on ${B2} with --broad option"
echo -e "\t macs2 callpeak -t ${bamfold}/${B2} -g hs -f BAMPE --broad -n ${P1} --outdir ${pkfold} -q 0.1 --keep-dup all 2>${pkfold}/${P1}.macs2Peak_summary.txt"

macs2 callpeak -t ${bamfold}/${B2} -g hs -f BAMPE --broad -n ${P1} --outdir ${pkfold} -q 0.1 --keep-dup all 2>${pkfold}/${P1}.macs2Peak_summary.txt

else

cd ${bgfold}
module load R/4.3.0

S1=$(echo $BG1 | sed 's/.bedgraph/.seacr.peaks/g')

echo "6. Calling SEACR Peaks on ${BG1}"
echo -e "\t bash $seacr ${bgfold}/${BG1} 0.01 non relaxed ${pkfold}/${S1}"

bash $seacr ${bgfold}/${BG1} 0.01 non relaxed ${pkfold}/${S1}

fi

echo -e "Peak calling complete \n"
###

##--------------------------------------------------------------## 9. Compute Heatmap

module load deeptools

echo "7. Computing heatmap of peaks for ${R1}"


if [ -z "${seacr}" ]
then

echo -e "\t using MACS2 peaks"
matfile=$(echo $BG1 | sed 's/.bedgraph/_MACS2.mat.gz/g')
M1=$(echo $BG1 | sed 's/.bedgraph/.macs2_summits.bed/g')

else

echo -e "\t using SEACR peaks"
M1=$(echo $S1 | sed 's/.seacr.peaks.relaxed.bed/.seacr.peaks.summitRegion.bed/g')
matfile=$(echo $S1 | sed 's/.seacr.peaks.relaxed.bed/_SEACR.mat.gz/g')
awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' ${pkfold}/${S1} > ${pkfold}/${M1}

fi

echo -e "\t computeMatrix reference-point -S ${covfold}/${BW2} -R ${pkfold}/${M1} --skipZeros -o ${pkfold}/${matfile} -p $threads -a 3000 -b 3000 --referencePoint center"
computeMatrix reference-point -S ${covfold}/${BW2} -R ${pkfold}/${M1} --skipZeros -o ${pkfold}/${matfile} -p $threads -a 3000 -b 3000 --referencePoint center

echo -e "\t plotHeatmap -m ${pkfold}/${matfile} -out ${rsfold}${matfile}_SEACR_heatmap.png --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${M1}""
plotHeatmap -m ${pkfold}/${matfile} -out ${rsfold}/${matfile}_heatmap.png --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${M1}"

echo -e " ******** All process have been completed ******** \n \n"
###