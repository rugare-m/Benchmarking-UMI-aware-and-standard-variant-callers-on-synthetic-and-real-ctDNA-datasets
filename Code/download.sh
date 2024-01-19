#!/bin/bash

#SBATCH --job-name=dietz_download           # Job name
#SBATCH --mem=32gb                   # Job memory request
#SBATCH --time=25:05:00              # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log     # Standard output and error log
pwd; hostname; date


declare -a arr=("SRR3401407" "SRR3401408" "SRR3401415" "SRR3401416" "SRR3401417" "SRR3401418")
#sratoolkit_dir="/users/rugarem/volatile/sratoolkit.3.0.0-centos_linux64/bin/"
fastq_dir="/users/rugarem/volatile/ctDNA/"

for srr in "${arr[@]}"
do 
    # prefetch to download the SRA data
    prefetch -p "$srr"
    
    # fasterq-dump to convert SRA to fastq
    fasterq-dump --include-technical -p -S "$srr"
    
    # Compress the resulting fastq files using gzip
    bgzip ./*.fastq
done

