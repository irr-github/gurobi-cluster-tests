#! /bin/bash    
#SBATCH --job-name=Snakemake
#SBATCH --output=slurmlogs/test-%j.out
#SBATCH --error=slurmlogs/test-%j.er
#SBATCH --mem=2048
#SBATCH --nodes=1-2
#SBATCH --cores=2
#SBATCH --qos=priority

PORT=1080
# sh workflow/scripts/kill-tunnel.sh $PORT
source /p/system/modulefiles/defaults/piam/module_load_piam
module load anaconda/2023.09
echo "Starting license server tunnel"
sh workflow/scripts/start-license-tunnel.sh $PORT
gurobi_cl --license &> "gurobi.log"
# snakemake --debug --cores 1 --use-conda
sh workflow/scripts/kill-tunnel.sh $PORT