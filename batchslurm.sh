#! /bin/bash
#SBATCH --ntasks 1     
#SBATCH --job-name=Snakemake
#SBATCH --output=slurmlogs/test-%j.out
#SBATCH --error=slurmlogs/test-%j.er
#SBATCH --mem=2048
#SBATCH --nodes=1 
#SBATCH --qos=io

snakemake --debug