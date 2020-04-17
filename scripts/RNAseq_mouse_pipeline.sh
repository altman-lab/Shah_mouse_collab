#!/bin/bash

##### Data cleaning RNA-seq data ######
########################################
## Setup
########################################
lsblk

sudo mkfs -t ext4 /dev/nvme1n1
sudo mkdir -p ~/genomics
sudo mount /dev/nvme1n1 ~/genomics
### Change permissions
sudo chmod 777 -R ~/genomics

########################################
## Upload data
########################################
## From local to S3
aws s3 sync /Users/kim/Documents/_Altman/Shah_collab/AM_mouse_chimera_TB/data_raw/fastq/ \
    s3://kadm-data/shah

## From S3 to EC2 instance    
## Using fuse
mkdir ~/genomics/data
sudo chmod 777 -R ~/genomics/data

s3fs kadm-data ~/genomics/data -o passwd_file=~/.passwd-s3fs \
    -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007
    ## Get uid and gid with `id`

mkdir ~/genomics/results
sudo chmod 777 -R ~/genomics/results

########################################
## Merge file from each sample
########################################
mkdir -p ~/genomics/results/fastq_merge
cd ~/genomics/data/shah

cd Sample_Lu1
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu01_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu01_R2.fastq.gz 
cd ..
cd Sample_Lu2
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu02_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu02_R2.fastq.gz 
cd ..
cd Sample_Lu3
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu03_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu03_R2.fastq.gz 
cd ..
cd Sample_Lu4
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu04_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu04_R2.fastq.gz 
cd ..
cd Sample_Lu5
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu05_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu05_R2.fastq.gz 
cd ..
cd Sample_Lu6
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu06_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu06_R2.fastq.gz 
cd ..
cd Sample_Lu7
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu07_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu07_R2.fastq.gz 
cd ..
cd Sample_Lu8
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu08_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu08_R2.fastq.gz 
cd ..
cd Sample_Lu9
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu09_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu09_R2.fastq.gz 
cd ..
cd Sample_Lu10
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu10_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu10_R2.fastq.gz 
cd ..
cd Sample_Lu11
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu11_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu11_R2.fastq.gz 
cd ..
cd Sample_Lu12
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu12_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu12_R2.fastq.gz 
cd ..
cd Sample_Lu13
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu13_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu13_R2.fastq.gz 
cd ..
cd Sample_Lu14
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu14_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu14_R2.fastq.gz 
cd ..
cd Sample_Lu15
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu15_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu15_R2.fastq.gz 
cd ..
cd Sample_Lu16
cat *R1_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu16_R1.fastq.gz 
cat *R2_00[0-9].fastq.gz > ~/genomics/results/fastq_merge/Lu16_R2.fastq.gz 
cd ..

cd 

########################################
## Quality assessment 1
## genomics/results/results_fastqc/fastqc_raw
########################################
## Assess read quality using FastQC

mkdir -p ~/genomics/results/results_fastqc/fastqc_raw/

for filename in ~/genomics/results/fastq_merge/*fastq.gz ;
do
echo "Starting FastQC analysis of" $filename

fastqc $filename \
       -o ~/genomics/results/results_fastqc/fastqc_raw/ \
       -t 15
done

## Download
aws s3 sync ~/genomics/results/ s3://shah-results/
## aws s3 sync s3://shah-results/results_fastqc ./results/results_fastqc
    
########################################
## Adapter removal
## genomics/data/fastq_trim
########################################
## Remove adapters 
## Remove reads with > 1 ambiguous base
## Trim ends until reach base with quality 30+
## Remove reads < 15 bp
mkdir -p ~/genomics/results/fastq_trim

paste <(ls ~/genomics/results/fastq_merge/*R1.fastq.gz) \
      <(ls ~/genomics/results/fastq_merge/*R2.fastq.gz) |

while read file1 file2;
do
  name1=$(paste -d '\0' \
            <(echo 'genomics/results/fastq_trim/') \
            <(awk -F'[_]' '{print $1}' <(basename $file1)))
  
  AdapterRemoval \
    --file1 $file1 \
    --file2 $file2 \
    --basename $name1 --gzip \
    --trim5p 15 --maxns 1 --minlength 15 \
    --trimqualities --minquality 30 \
    --threads 15

done
    
########################################
## Quality assessment 2
## genomics/results/results_fastqc/fastqc_trim
########################################
## Assess trimmed read quality using FastQC.

mkdir -p ~/genomics/results/results_fastqc/fastqc_trim/

for filename in ~/genomics/results/fastq_trim/*pair[12].truncated.gz ;
do
echo "Starting FastQC analysis of" $filename
fastqc $filename \
       -o ~/genomics/results/results_fastqc/fastqc_trim/ \
       -t 15
done
    
## Download
aws s3 sync ~/genomics/results/ s3://shah-results/
#aws s3 sync s3://shah-results/ ./results

########################################
## Alignment
## genomics/data/bam
########################################
## Get ref data files
mkdir -p ~/genomics/results/STARindex.mouse
mkdir -p ~/genomics/results/STARref.mouse
cd ~/genomics/results/STARref.mouse

sudo curl -O ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
gunzip Mus_musculus.GRCm38.99.gtf.gz

sudo curl -O ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

cd

## Make genome index
STAR --runMode genomeGenerate \
     --genomeDir ~/genomics/results/STARindex.mouse \
     --genomeFastaFiles ~/genomics/results/STARref.mouse/Mus_musculus.GRCm38.dna.primary_assembly.fa \
     --sjdbGTFfile ~/genomics/results/STARref.mouse/Mus_musculus.GRCm38.99.gtf \
     --sjdbOverhang 99 \
     --runThreadN 15
     
## Align with STAR
mkdir -p ~/genomics/results/bam/

paste <(ls ~/genomics/results/fastq_trim/*pair1.truncated.gz) \
      <(ls ~/genomics/results/fastq_trim/*pair2.truncated.gz) |

while read file1 file2;
do
    echo "Aligning" $(basename  -- "$file1");
    
    name=$(paste -d '\0' \
            <(echo 'genomics/results/bam/') \
            <(awk -F'[_.]pair' '{print $1}' <(basename $file1)) \
            <(echo '_'))

    STAR --genomeDir genomics/results/STARindex.mouse \
         --readFilesIn $file1 $file2 \
         --readFilesCommand zcat \
         --outFileNamePrefix $name \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 15 \
         --runRNGseed 8756
done

########################################
## Assess aligned reads
## genomics/results/results_cleaning
########################################
## Get refFlat genome for Picard
mkdir -p ~/genomics/results/PICARDref
cd ~/genomics/results/PICARDref
sudo curl -O http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/refFlat.txt.gz
gunzip -c refFlat.txt.gz > refFlat.mouse.txt
## Remove chr in chromosome name to match ensembl alignment
sed 's/chr//' refFlat.mouse.txt > refFlat.mouse.ensembl.txt
cd

## median CV of gene model coverage
mkdir -p ~/genomics/results/results_cleaning/

for bam_file in ~/genomics/results/bam/*sortedByCoord.out.bam ;
do
    java -XX:ParallelGCThreads=15 \
        -jar ~/genomics/anaconda3/share/picard-2.22.1-0/picard.jar \
        CollectRnaSeqMetrics \
        REF_FLAT=~/genomics/results/PICARDref/refFlat.mouse.ensembl.txt \
        INPUT=$bam_file  \
        OUTPUT=~/genomics/results/results_cleaning/temp.tsv \
        ASSUME_SORTED=true \
        STRAND_SPECIFICITY=NONE \
        MINIMUM_LENGTH=500

    #Append results
    echo $bam_file >> ~/genomics/results/results_cleaning/bam.metrics.tsv
    cat ~/genomics/results/results_cleaning/temp.tsv >> ~/genomics/results/results_cleaning/bam.metrics.tsv
    #Remove this iteration
    rm ~/genomics/results/results_cleaning/temp.tsv
done
        
## mapped_reads_w_dups aka alignment percentage

for bam_file in ~/genomics/results/bam/*sortedByCoord.out.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> ~/genomics/results/results_cleaning/summary.alignment.tsv
    
    samtools flagstat -@ 10 $bam_file \
    >> ~/genomics/results/results_cleaning/summary.alignment.tsv
done

aws s3 sync ~/genomics/results/ s3://shah-results/
## aws s3 sync s3://shah-results/results_cleaning ./results/results_cleaning

########################################
## Quality filter BAM
## genomics/data/bam_filter
########################################
## -h: keep header
## -f 3: keeps paired reads where both mapped
## -F 1284: removes unmapped reads, non-primary alignments, and PCR duplicates
## -q 30: min MAPQ of 30

## Paired reads
mkdir -p genomics/results/bam_filter_paired/

for bam_file in genomics/results/bam/*sortedByCoord.out.bam ;
do
  filename1=$(paste -d '\0' \
            <(echo 'genomics/results/bam_filter_paired/') \
            <(awk -F'[_]Aligned' '{print $1}' <(basename $bam_file)) \
            <(echo '_filter_paired.bam'))
  
  echo "Filtering" $bam_file          
  samtools view $bam_file \
      -h -f 3 -F 1284 -q 30 \
      -@ 15 \
      > $filename1
done


## alignment percentage

for bam_file in genomics/results/bam_filter_paired/*filter_paired.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> genomics/results/results_cleaning/summary.align.filter.paired.tsv
    
    samtools flagstat -@ 15 $bam_file \
    >> genomics/results/results_cleaning/summary.align.filter.paired.tsv
done

########################################
## Count reads in genes
## genomics/data/counts
########################################
## List all possible features to count in annotation file
#cut -f 3 genomics/results/STARref.mouse/Mus_musculus.GRCm38.99.gtf | sort | uniq
mkdir -p genomics/results/counts

## Count reads in genes
featureCounts -T 14 -g gene_id -t exon -p \
  -a genomics/results/STARref.mouse/Mus_musculus.GRCm38.99.gtf \
  -o genomics/results/counts/Shah.featurecounts.paired.tsv \
  genomics/results/bam_filter_paired/*filter_paired.bam

########################################
## Save to S3 storage
########################################

aws s3 sync ~/genomics/results/ s3://shah-results
#aws s3 sync s3://shah-results/results_cleaning ./results/results_cleaning
#aws s3 sync s3://shah-results/counts ./results/counts
#aws s3 sync s3://shah-results ./results

################# END ##################