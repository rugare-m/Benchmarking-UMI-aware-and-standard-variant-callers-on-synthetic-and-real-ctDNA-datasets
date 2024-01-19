#!/bin/bash

# Job name:
#SBATCH --job-name=callers
#SBATCH --mem=64G
# Wall clock limit:
#SBATCH --time=24:00:30

time python bcftools.py
time python lofreq.py
time python freebayes.py
time python mutect2.py
