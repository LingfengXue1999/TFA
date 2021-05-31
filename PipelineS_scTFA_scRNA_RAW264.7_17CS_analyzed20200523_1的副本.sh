#!/bin/bash
source /home/yihanlin_pkuhpc/lustre2/wuyan/.bashrc; 
cd /home/yihanlin_pkuhpc/lustre3/wuyan/scTFA/raw264.7_17CS/data;

n_threads=10;

# build index for mapping
mkdir STAR_index ;
STAR  --runMode genomeGenerate --genomeDir ./STAR_index --runThreadN ${n_threads} \
--genomeFastaFiles /home/yihanlin_pkuhpc/lustre2/wuyan/resource/mousegenome/GRCm38.p6.genome.fa \
--sjdbGTFfile /home/yihanlin_pkuhpc/lustre2/wuyan/resource/mousegenome/gencode.vM25.annotation.gtf --sjdbOverhang 75 ;

for data in `ls | grep "SRR"`
do
  {
  cd ${data};
  # transform sra files to fastq files
  fastq-dump --split-files ${data}.sra;

  mkdir ./map;

  # trim with Trimmomatic
  java -jar /home/yihanlin_pkuhpc/lustre2/wuyan/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads ${n_threads} -trimlog ./map/trim.log \
  ${data}_1.fastq ${data}_2.fastq \
  ./map/${data}_output_forward_paired_R1.fq.gz ./map/${data}_output_forward_unpaired_R1.fq.gz \
  ./map/${data}_output_reverse_paired_R2.fq.gz ./map/${data}_output_reverse_unpaired_R2.fq.gz \
  ILLUMINACLIP:/home/yihanlin_pkuhpc/lustre2/wuyan/app/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:36 ;

  # align with STAR
  STAR \
      --runMode alignReads \
      --outFileNamePrefix ./map/STAR \
      --genomeDir ../STAR_index \
      --runThreadN ${n_threads} \
      --readFilesIn ./map/${data}_output_forward_paired_R1.fq.gz ./map/${data}_output_reverse_paired_R2.fq.gz \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outBAMsortingThreadN ${n_threads} ;

  # count with TPMCalculator
  TPMCalculator -g /home/yihanlin_pkuhpc/lustre2/wuyan/resource/mousegenome/gencode.vM25.annotation.gtf \
  -b /home/yihanlin_pkuhpc/lustre3/wuyan/scTFA/raw264.7_17CS/data/velocyto/${data}/${data}.bam -c 75 -p ;

  cd ..;
} &
wait
done
