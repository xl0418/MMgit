#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=dataanalysis
#SBATCH --output=job-%j.log
#SBATCH --mem=1000
module load R/3.2.3-foss-2016a
time Rscript test2-1t.R $1
