# PIK Cluster tests

Simple tests to see whether we cna get gurobi to work on the PIK cluster with snakemake

# Setup

Gurobi license activation from the compute nodes requries internet access. Workaround is an ssh tunnel to the compute nodes, which can be set-up on the compute nodes with
```
# interactive session on the compute nodes
srun --qos=priority --pty bash
# key pair gen (here ed25518 but can be rsa)
ssh-keygen -t ed25519 -f ~/.ssh/id_rsa.cluster_internal_exchange -C "$USER@cluster_internal_exchange"
# leave the compute nodes
exit
```
You will then need to add the contents of the public key `~/.ssh/id_rsa.cluster_internal_exchange.pub` to your authorised `~/.ssh/authorized_keys`

In addition you should have your .profile setup as per https://gitlab.pik-potsdam.de/rse/rsewiki/-/wikis/Cluster-Access
and add `module load anaconda/2024.10` (or latest) to it 

# General usage
- the run_snakemake.sh is the main. It should be run as `sbatch run_snakemake.sh` from the login nodes of the cluster.
    - snakemake takes care of further sbatch submissions but the gurobi license tunnel needs to be run from the compute nodes
- the Snakefile configures the execution of the snakemake workflow
- the config/pik_hpc_profile/config.yaml takes care of the ressource &qos allocation for each rule

# issues:
- rules not running fully as expected
  
## Challenges

> **License**: 
> - proper options need to be passed in order for the WSL license to be recognised. Otherwise you get a site mismatch id
> *Best solution is to set the GUROBI env parameters, otherwise you may need to explicitly pass the parameters to an env
> - license checking requires a connection to the internet, which is not available on the compute nodes
> - number of cores allowed on the license is limited
> - number of simultaneous gurobi instances is limited
