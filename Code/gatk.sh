#!/bin/bash

# Job name:
#SBATCH --job-name=gatk4
#SBATCH --mem=64G       
#SBATCH --time=24:00:30

python gatk.py