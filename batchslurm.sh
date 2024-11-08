#! /bin/bash
#SBATCH --ntasks 1     
#SBATCH --job-name=Snakemake
#SBATCH --output=slurmlogs/test-%j.out
#SBATCH --error=slurmlogs/test-%j.er
#SBATCH --mem=2048
#SBATCH --nodes=1-2
#SBATCH --cores=2
#SBATCH --qos=priority

PORT=1080
source /p/system/modulefiles/defaults/piam/module_load_piam
module load anaconda/2023.09
sh workflow/scripts/start-license-tunnel.sh $PORT
snakemake --debug --cores 1 --use-conda
sh workflow/kill-tunnel.sh $PORT