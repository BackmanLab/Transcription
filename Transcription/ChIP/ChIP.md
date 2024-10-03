## ChIP-seq reanalysis 

PMID32616013_ChIP is main dir for project. This project is reananysis of ChIP data for Pol 1, 2, and 3 in mouse ESCs. Find more info [here](https://pubmed.ncbi.nlm.nih.gov/32616013/)

all SRRs can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145874)

### Get FASTQs

Sample header:
```shell
#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=16 
#SBATCH --mem=100gb
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --time 48:00:00
#SBATCH --job-name=Samtools
#SBATCH --output=outlog
#SBATCH --error=errlog

```

Get FASTQS for ChIP data. [Fastq-dump tutorial](https://rnnh.github.io/bioinfo-notebook/docs/fastq-dump.html)

Download the following SRA files:
                        Run	                    # of Spots	# of Bases	Size	Published
- WT_Pol1_ChIP_B1T2 <kbd>SRR11150185</kbd>		49,069,244	14.7G	6Gb	        2020-06-03
- WT_Pol1_ChIP_B1T1 <kbd>SRR11150184</kbd>		48,744,574	14.6G	6.1Gb	    2020-06-03
- WT_Pol1_ChIP_B2T1 <kbd>SRR11150186</kbd>		35,624,523	10.7G	4.4Gb	    2020-06-03
- WT_Pol1_ChIP_B2T2 <kbd>SRR11150187</kbd>		32,688,512	9.8G	4.1Gb	    2020-06-03
- WT_Pol2_ChIP_B1T1 <kbd>SRR11150188</kbd>		29,659,368	8.9G	3.7Gb	    2020-06-03
- WT_Pol2_ChIP_B1T2 <kbd>SRR11150189</kbd>	    43,241,389	13G	    5.4Gb	    2020-06-03
- WT_Pol2_ChIP_B2T1 <kbd>SRR11150190</kbd>      42,011,385	12.6G	5.5Gb	    2020-06-03
- WT_Pol2_ChIP_B2T2 <kbd>SRR11150191</kbd>      46,504,618	14G	    6Gb	        2020-06-03
- WT_Input_rep1     <kbd>SRR11150182</kbd> 	    36,124,283	10.8G	4.5Gb	    2020-06-03
- WT_Input_rep2     <kbd>SRR11150183</kbd>      37,199,087	11.2G	4.7Gb	    2020-06-03

We will need the SRA toolkit to download the data. SRA is a data compression format for FASTQ sequencing data that makes data more portable and reduces it's memory footprint. Load the <kbd>sratoolkit module<kbd>

```shell

# Set working directory
cd /projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP

module load sratoolkit
```

```shell
fastq-dump --gzip --split-3 SRR11150185 SRR11150184 SRR11150186 SRR11150187 SRR11150188 SRR11150189 SRR11150190 SRR11150191 SRR11150182 SRR11150183
```
Faster version: 

```shell
##--SLURM header here--## 

# Load necessary modules
module load sratoolkit

# Set working directory
cd /projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP/fastq

# Define array of SRR numbers
SRRfiles=(SRR11150185 SRR11150184 SRR11150186 SRR11150187 SRR11150188 SRR11150189 SRR11150190 SRR11150191 SRR11150182 SRR11150183)

## define outfolder variable where FASTQs will go
outfold="$(dirname "$(pwd)")/${PWD##*/}"
echo "Fetching data to $outfold "

# For loop retrieves SRR FASTQs and change names for compatibility with
for SRR in ${SRRfiles[@]}; do
    echo "Fetching FASTQs from $SRR"
    fasterq-dump -O $outfold -e 16 --split-3 $SRR
    echo "$SRR retrieved"
done

# GZIP all FASTQs to save space
for file in *fastq; do
    gzip $file;
    echo "$file has been compressed"
done

# For loop retrieves SRR FASTQs and change names for compatibility with
for SRR in ${SRRfiles[@]}; do
      echo "Formatting $SRR FASTQ"

      for i in {1..2}; do
        mv ${SRR}_${i}.fastq.gz ${SRR}_R${i}.fastq.gz;
     done
done

```

Reference genome is built in this paper: [Construction and validation of customized genomes for human and mouse ribosomal DNA mapping](https://www.jbc.org/article/S0021-9258(23)01794-5/fulltext)
Get reference genome that contains rDNA repeats from [here](https://github.com/vikramparalkar/rDNA-Mapping-Genomes)

```shell
wget https://github.com/vikramparalkar/rDNA-Mapping-Genomes/raw/main/Mouse_mm39-rDNA_genome_v1.0.tar.gz?download=
```

hg38-rDNA_v1.0.fa /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/refR/
full Slurm commands can be found [here](https://slurm.schedmd.com/scancel.html). Wrote and ran the following FASTQC slurm submission script to assess pre-trimming quality of reads:

```shell
##--SLURM header here--## 

# load modules
module load fastqc

# The command to execute:
fastqc -t 16 *.fastq.gz

#Set your working directory
cd /projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP/fastq

mkdir Fastqc_Results

mv *.html Fastqc_Results
mv *fastqc.zip Fastqc_Results
```

### Bowtie Indices

Build Bowtie indices for reference genome prior to alignment  <kbd> bowtie2-build path_to_reference_genome.fa prefix_to_name_indexes </kbd> using default bowtie2 on quest, bowtie2/2.4.5

```shell
##--SLURM header here--## 

cd /projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP/genome/Mouse_mm39-rDNA_genome_v1.0

module load bowtie2/2.4.5

bowtie2-build /projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP/genome/Mouse_mm39-rDNA_genome_v1.0/mm39.fa mm39
```

After running <kbd> bowtie_build.sh </kbd>, you should have the following files in your index directory:
```
mm39.1.bt2  mm39.2.bt2  mm39.3.bt2  mm39.4.bt2  mm39.fa  mm39.fa.fai  mm39.rev.1.bt2  mm39.rev.2.bt2
```

### Next, we run Trimmomatic on the the FASTQS

Reads needed to be trimmed to remove adapter sequences from Nextera/Illumina sequencing and low quality bases. [Trimmomatic ](http://www.usadellab.org/cms/?page=trimmomatic)for this purpose. Run FASTQC first to get the overrepresented sequences to get the adapters that were used and need to be clipped. The HTMLs will have these sequences in them. using FASTQ, i found that the adapters in my reads matched TruSeq3-PE-2.fa adapter file provided by Trimmomatic. Adapter sequences are [here](https://github.com/timflutre/trimmomatic/blob/master/adapters/). Used the following tutorials to guide scripting and parameters:

1. [data carpentry](https://datacarpentry.org/wrangling-genomics/03-trimming.html)
2. [Trimmomatic Manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
3. [A data scientist's blog](https://sites.psu.edu/biomonika/2015/07/03/trim-now-or-never-pe-and-mp-adaptor-removal/)
4. These Biostars entries: [1](https://www.biostars.org/p/257586/) & [2](https://www.biostars.org/p/366041/)

**Wrote the following Slurm submission script for trimming:**

```shell
##--SLURM header here--## 

#Set your working directory
cd /projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP/fastq

#load modules
module load trimmomatic/0.39

#declare function
function trim_data {
R1=$1;
R2=$(echo $R1 | sed 's/_R1/_R2/g')
Out1=$(echo $R1 | sed 's/.fastq.gz/./g')
Out2=$(echo $R2 | sed 's/.fastq.gz/./g')

java -jar /software/trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 -trimlog ${Out1}_trimlog \
$R1 $R2 ${Out1}trimmed.fastq.gz ${Out1}Fsingle.fastq.gz ${Out2}trimmed.fastq.gz ${Out2}Rsingle.fastq.gz \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30;

rm ${Out1}Fsingle.fastq ${Out2}Rsingle.fastq;
}

for file in *R1.fastq.gz;
do trim_data $file;
done

mkdir Trimmed_Reads

mv *trimmed.fastq.gz Trimmed_Reads
mv *_trimlog Trimmed_Reads

```

### Map reads

Check these two HBC references for more information on mapping reads using Bowtie2

1. [the first one](https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/04_alignment_using_bowtie2.html)
2. [the second one](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html)

Here is the mapping script for trimmed reads

```shell
##--SLURM header here--## 

# Load modules
module load bowtie2/2.4.5

# Set working directory
cd /projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP/fastq/Trimmed_Reads

echo "$(date): Processing alignment of FASTQ sequence files"

echo "Generating directories."

genomefold=/projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP/genome/Mouse_mm39-rDNA_genome_v1.0/mm39

echo "Bowtie genome assembly is located in $genomefold"

outfold="$(dirname "$(pwd)")/${PWD##*/}/output"
[ ! -d $outfold ]&&mkdir $outfold

# Define function align_data
function align_data {

R1=$1;
R2=$(echo $R1 | sed 's/_R1/_R2/g')
Out=$(echo $R1 | sed 's/_R1.trimmed.fastq.gz/./g')

bowtie2 -p 16 -q -X 2000 -x $genomefold -1 ${R1} -2 ${R2} -S "$outfold/${Out}sam"


}

echo "Aligning files"
for file in *R1.trimmed.fastq.gz;
do echo "${file}";
align_data $file;
done

```

### SAM to BAM

Here's the SAM to BAM and coverage script, based on suggestions from the [rDNA alignment paper](https://www.jbc.org/article/S0021-9258(23)01794-5/fulltext)

```shell

##--SLURM header here--## 

cd /projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP/fastq/Trimmed_Reads/output

module load samtools/1.6
module load deeptools/3.1.1

# Make QC directory
qcfold="$(dirname "$(pwd)")/${PWD##*/}_qc"
[ ! -d $qcfold ]&&mkdir $qcfold

# Make coverage directory
covfold="$(dirname "$(pwd)")/${PWD##*/}_coverage"
[ ! -d $covfold ]&&mkdir $covfold

# Make .BAM directory
bamfold="$(dirname "$(pwd)")/${PWD##*/}_BAMs"
[ ! -d $bamfold ]&&mkdir $bamfold

# Define function sam_tools
function sam_tools {
SAM=$1;
BAM=$(echo $SAM | sed 's/.sam/.bam/g')
BAMsort=$(echo $BAM | sed 's/.bam/.sort.bam/g')
BW=$(echo $BAMsort | sed 's/.sort.bam/./g')

samtools view -h -F 4 -q 1 -bS $SAM > $BAM
samtools sort $BAM -o $BAMsort
samtools index $BAMsort


# Bam coverage to generate BigWig
bamCoverage --bam $BAMsort --normalizeUsing CPM --outFileName "${BW}-.bigWig" --filterRNAstrand reverse --binSize 1 --numberOfProcessors 60
bamCoverage --bam $BAMsort --normalizeUsing CPM --outFileName "${BW}+.bigWig" --filterRNAstrand forward --binSize 1 --numberOfProcessors 60

}

echo "Sorting alignments and generating BigWigs"
for file in *.sam;
do echo "${file}";
sam_tools $file;
done

echo "Gathering statistics on BAM files"
for file in *sort.bam;
do echo "$file";

samtools flagstat $file > \
${file}.FLAGSTAT.txt

done
echo "Moving processed files to respective directories."

mv *FLAGSTAT.txt $qcfold
mv *.bam $bamfold
mv *.bigWig $covfold

echo "FASTQ processing complete!"


```

### Call Peaks

Next, we call peaks on the BAMs using [MACS2](https://pypi.org/project/MACS2/). The methods from  [Genome-wide analyses of chromatin interactions after the loss of Pol I, Pol II, and Pol III ](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02067-3#Sec13) state that: 

"Pooled ChIP libraries were prepared and sequenced at a sequencing depth of ~ 90–100 million reads per sample. Raw fastq reads were trimmed by trim_galore and mapped to the mouse reference genome mm10 using bowtie2 (version 2.3.2) [85] with the default parameters “best –k 1 –m 1 and –l 18.” Unmapped reads, low-quality mapped reads, and PCR duplicates were discarded. Only uniquely mapped data were retained for the downstream analysis. In the next step, we carried out peak calling individually for each replicate against the input control using MACS2 (version 2.1.0) with the “-c” option and a p value threshold of 10−5 to ensure high confidence [86]. The peaks that overlapped a peak from the other replicate of the same RNAP sample by at least 1 bp were retained. In these cases, the new peak equaled the combined coordinates of all the overlapping peaks, considering both the number of biological replicates and the treatment."

From their GEO deposit:

"Peaks were called against input control using the MACS2 function callpeak on each replicate with options ‘-g mm -q 0.0001 --nomodel --extsize 200’ and p-value threshold of 10-5 to ensure high confidence. Common peaks between the replicates were used for downstream analysis."

I am following the [MACS2 tutorial here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html). They recommend the following: 

-t: The input data file (this is the only REQUIRED parameter for MACS). This is the treatment (the sorted BAM)
-c: The control or mock data file (WT HCT116 input controls: SRR11150182 and SRR11150183). These files are non ChIPped genomic background

From the github it looks like we can pool the two input controls and they will be merged or pooled

-t/--treatment FILENAME

This is the only REQUIRED parameter for MACS. The file can be in any supported format -- see detail in the --format option. If you have more than one alignment file, you can specify them as -t A B C. MACS will pool up all these files together.
-c/--control

The control, genomic input or mock IP data file. Please follow the same direction as for -t/--treatment. 

Running into an error with MACS2 indexing. Going to try merging input BAMs and see if that helps:

```shell

#Merge input BAMS
# echo "Merging input BAM files for pseudoreplicates..."
# samtools merge -u ${tmpDir}/${NAME1}input_${NAME2}input_merged.bam $baseDir/${inputFile1} $baseDir/${inputFile2}

#Merge input BAMS
samtools merge -u input_merged.bam SRR11150182.sort.bam SRR11150183.sort.bam

```
I used a less stringent cutoff - an FDR (q value) of 0.01, for calling peaks

```shell

##--SLURM header here--## 

cd /projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP/fastq/Trimmed_Reads/output_BAMs

module load MACS2/2.1.0

# Make output directory
peaksfold="$(dirname "$(pwd)")/peaks"
[ ! -d $peaksfold ]&&mkdir $peaksfold

function callPeaks {

## Define variables
C=input_merged.bam # input control r
sortedBAM=$1; # sorted bam
Out=$(echo $sortedBAM | sed 's/.sort.bam/./g') # output root

echo "processing $sortedBAM"

## Call MACS
macs2 callpeak -t $sortedBAM -c $C -f BAMPE -n "$peaksfold/${Out}" -g mm --nomodel -q 1e-2 --keep-dup all --extsize 200 # can use -q 1e-2 or -p 1e-5 if q values are too stringent

}

function sortMACS {

sortedBAM=$1
base=`basename $sortedBAM .sort.bam`

 mkdir $base
 mv $(dirname "$(pwd)")/${base}.narrowPeak $peaksfold
 mv $(dirname "$(pwd)")/${base}.bed $peaksfold
 mv $(dirname "$(pwd)")/${base}.xls $peaksfold

}

for i in *.sort.bam;
do callPeaks $i;
done

for i in *.sort.bam;
do sortMACS $i;
done
```
Checked number of peaks called for background files and there only a couple hundred, validating this as a control.


```
[lmc0633@quser33 peaks]$ wc -l SRR11150182._peaks.narrowPeak
218 SRR11150182._peaks.narrowPeak
[lmc0633@quser33 peaks]$ wc -l SRR11150183._peaks.narrowPeak
244 SRR11150183._peaks.narrowPeak
[lmc0633@quser33 peaks]$ 
```

I now want to merge my peak files. There are a few different ways of merging technical and biological replicates. One is to use [Bedtools Intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) to build a union of peaks. This is described for ChIP-seq data [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/handling-replicates-bedtools.html). A more comprehensive quality control approach to merging replicates is described [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html). Since this data is already validated elsewhere and we're just remapping, I'm going to proceed without IDR. Here is the Bedtools intersect approach.

```shell

module load bedtools

## First, concatenate the two technical replicates for each biological replicate into a single bed file and sort them by chromosome and position.
cat SRR11150184._peaks.narrowPeak SRR11150185._peaks.narrowPeak | sort -k1,1 -k2,2n  > pol1.t1.narrowPeak

cat SRR11150186._peaks.narrowPeak SRR11150187._peaks.narrowPeak | sort -k1,1 -k2,2n  > pol1.t2.narrowPeak

## intersect these if they overlap by %5 (-f 0.05) and Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported. Restricted by -f and -r (-wo)
bedtools intersect -f 0.05 -r -a pol1.t1.narrowPeak -b pol1.t2.narrowPeak -wo > pol1.union.bed

```

Alternatively, there is a method proposed [here]https://ro-che.info/articles/2018-07-11-chip-seq-consensus) that uses genomic ranges in R to find peaks overlapping in more than one replicate (2) and removes peaks from the dataset that only occur in a single replicate. This seems more robust. I will use this

#### Later modification

I looked into merging replicates using Coverage() in R and unfortunately, you lose all the metadate (qValue, pValue, signal) when doing so. I decided to merge my technical replicates at the BAM level and then call the peaks again on those. Following calling peaks, I plan to use Bedtools intersect to merge the biological replicates, since it preserves the metadata and other information

```shell
##--SLURM header here--## 

cd /projects/b1042/BackmanLab/Lucas/PMID32616013_ChIP/fastq/Trimmed_Reads/output_BAMs

module load samtools

#Merge input BAMS
samtools merge -u pol1.rep1.sort.bam SRR11150184.sort.bam SRR11150185.sort.bam

#Merge input BAMS
samtools merge -u pol1.rep2.sort.bam SRR11150186.sort.bam SRR11150187.sort.bam

#Merge input BAMS
samtools merge -u pol2.rep1.sort.bam SRR11150188.sort.bam SRR11150189.sort.bam

#Merge input BAMS
samtools merge -u pol2.rep2.sort.bam SRR11150190.sort.bam SRR11150191.sort.bam

```

```shell

module load bedtools

## intersect these if they overlap by %5 (-f 0.05) and Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported. Restricted by -f and -r (-wo)
bedtools intersect -r -wo -f 0.1 -a pol1.rep1._peaks.narrowPeak  -b pol1.rep2._peaks.narrowPeak  > pol1.bedtools.union.bed

bedtools intersect -r -wo -f 0.1 -a pol2.rep1._peaks.narrowPeak  -b pol2.rep2._peaks.narrowPeak  > pol2.bedtools.union.bed

```

After intersecting peaks, peak files and coverage files were analyzed in R. 