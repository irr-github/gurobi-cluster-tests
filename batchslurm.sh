#! /bin/bash    
#SBATCH --export=PATH,LD_LIBRARY_PATH,GUROBI_HOME
#SBATCH --job-name=Snakemake
#SBATCH --output=slurmlogs/test-%j.out
#SBATCH --error=slurmlogs/test-%j.er
#SBATCH --mem=2048
#SBATCH --qos=priority
#SBATCH --cpus-per-task=4 

#set up the tunnel to the login nodes (does not work well as a subprocess, even with source)
PORT=1080
ssh -fN -D $PORT $USER@login01 &

export https_proxy=socks5://127.0.0.1:$PORT
export SSL_CERT_FILE=/p/projects/rd3mod/ssl/ca-bundle.pem_2022-02-08
export GRB_CAFILE=/p/projects/rd3mod/ssl/ca-bundle.pem_2022-02-08

# Add gurobi references
export GUROBI_HOME="/p/projects/rd3mod/gurobi1103/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
export GRB_LICENSE_FILE=/p/projects/rd3mod/gurobi_rc/gurobi.lic
export GRB_CURLVERBOSE=1

# avoid GUROBIPY complaining
export NUMEXPR_MAX_THREADS=16

# # check the license (FAILS due to core count)
gurobi_cl --license &> logs/gurobi.log

# PROFILE the machine (what gurobi license check does under the hood)
# grbprobe &> gurobi.log

snakemake --use-conda --debug --cores 1 --forcerun log_test