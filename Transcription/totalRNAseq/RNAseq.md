
## Total RNAseq alignment scripts

Below are example scripts used to align and count the total RNA seq data generated for this paper

### Step 1: Build references for STAR:

Example header
```shell 
#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 48:00:00
#SBATCH --job-name=STAR index

```

```shell

##--SLURM header here--##

module load STAR/2.6.0

cd /projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/STAR_107

STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir STAR \
--genomeFastaFiles /projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/STAR_107/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/STAR_107/Homo_sapiens.GRCh38.107.gtf \
--sjdbOverhang 74

```

### Step 2: Build RSEM references:

```shell
##--SLURM header here--##

cd /projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/Gencode

# load modules you need to use
module load rsem/1.3.0
module load STAR/2.5.2

# A command you actually want to execute:
rsem-prepare-reference -p 16 --gtf gencode.v41.primary_assembly.annotation.gtf --star GRCh38.primary_assembly.genome.fa /projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/Gencode/RSEM

```

### Step 3: Align and count:

```shell

##--SLURM header here--##

# Load modules
module load python/anaconda3.6
module load STAR/2.6.0
module load rsem/1.3.3
module load samtools/1.6
module load deeptools/3.1.1
module load fastqc

# Set working directory
cd /projects/b1042/BackmanLab/Lucas/totalRNAseq/

echo "$(date): Processing alignment of FASTQ sequence files"

echo "Generating directories."
# make sure the STAR and RSEM references are index references
genomefold=/projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/Gencode/STAR
gtffold=/projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/Gencode/gencode.v41.annotation.gtf
RSEMfold=/projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/Gencode/RSEM/RSEM

echo "STAR genome assembly is located in $genomefold"
echo "Annotations file is located in $gtffold"
echo "RSEM reference files are located in $RSEMfold"

outfold="$(dirname "$(pwd)")/${PWD##*/}/output"
[ ! -d $outfold ]&&mkdir $outfold

# Make QC directory
qcfold="$(dirname "$(pwd)")/${PWD##*/}/fastqc"
[ ! -d $qcfold ]&&mkdir $qcfold

# Make coverage directory
sortfold="$(dirname "$(pwd)")/${PWD##*/}/coverage"
[ ! -d $sortfold ]&&mkdir $sortfold

# Make .BAM directory
bamfold="$(dirname "$(pwd)")/${PWD##*/}/BAMs"
[ ! -d $bamfold ]&&mkdir $bamfold

# Make .OUT directory
outfold2="$(dirname "$(pwd)")/${PWD##*/}/OUT"
[ ! -d $outfold2 ]&&mkdir $outfold2

# Make .HTSEQ_COUNTS directory
cntfold="$(dirname "$(pwd)")/${PWD##*/}_HTSEQ_counts"
[ ! -d $cntfold ]&&mkdir $cntfold

# Make .RSEM_ COUNTS directory
cntfold2="$(dirname "$(pwd)")/${PWD##*/}_RSEM_counts"
[ ! -d $cntfold2 ]&&mkdir $cntfold2

# Define function align_data
function align_data {

R1=$1;
R2=$(echo $R1 | sed 's/R1/R2/g')

STAR --runThreadN 24 --quantMode TranscriptomeSAM --genomeDir $genomefold --readFilesIn ${R1} ${R2} --readFilesCommand zcat --outFileNamePrefix "$outfold/${R1%.*.*}_";
}

# Define function sam_tools
function sam_tools {
R1=$1;
BAM=$(echo $R1 | sed 's/_R1_001_Aligned.out.sam/_Aligned.out.bam/g')
BAMsort=$(echo $BAM | sed 's/_Aligned.out.bam/_Aligned.out.sort.bam/g')
BW=$(echo $BAMsort | sed 's/_Aligned.out.sort.bam/./g')

samtools view -b $R1 > $BAM
samtools sort $BAM > $BAMsort
samtools index $BAMsort

# Bam coverage to generate BigWig
bamCoverage --bam $BAMsort --normalizeUsing CPM --outFileName "${BW}-.bigWig" --filterRNAstrand reverse --binSize 1 --numberOfProcessors 60
bamCoverage --bam $BAMsort --normalizeUsing CPM --outFileName "${BW}+.bigWig" --filterRNAstrand forward --binSize 1 --numberOfProcessors 60

}

echo "Aligning files"
for file in *R1_001.fastq.gz;
do echo "${file}";
align_data $file;
done

echo "Running Fastqc on all FASTA files"
fastqc -t 12 *.fastq.gz
mv *.html $qcfold
mv *fastqc.zip $qcfold

# Change directory to $Outfold
cd $outfold

echo "Sorting alignments and generating BigWigs"
for file in *Aligned.out.sam;
do echo "${file}";
sam_tools $file;
done

echo "Gathering statistics on BAM files"
for file in *sort.bam;
do echo "$file";

samtools flagstat $file > \
${file}.FLAGSTAT.txt

done

echo "Generating HTSEQ count files from sample alignments"
for file in *.sam;
do echo "${file}";
output_file=$(echo $file | sed 's/_R1_001_Aligned.out.sam/.counts/g')
htseq-count -f sam -r pos -s no -t exon -i gene_id -m union "${file}" $gtffold > "${output_file}"
done

echo "Generating RSEM count files from sample alignments"
for file in *.toTranscriptome.out.bam;
do echo "${file}";
output_file=$(echo $file | sed 's/_R1_001_Aligned.toTranscriptome.out.bam/.rsem/g')
rsem-calculate-expression -p 24 --alignments --paired-end "${file}" $RSEMfold "${output_file}"
done

echo "Moving processed files to respective directories."

mv *FLAGSTAT.txt $qcfold
mv *.bam $bamfold
mv *.out $outfold2
mv *.counts $cntfold
mv *.results $cntfold2
mv *.bigWig $sortfold

echo "FASTQ processing complete!"

```

