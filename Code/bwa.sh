#!/bin/bash

# Job name:
#SBATCH --job-name=bwa-mem2
#SBATCH --mem=128G           # total memory per node
# Wall clock limit:
#SBATCH --time=24:00:30


python bwa.py
