# PIK Cluster tests

Simple tests to see whether we cna get gurobi to work on the PIK cluster with snakemake

## Challenges

> **License**: 
> - proper options need to be passed in order for the WSL license to be recognised. Otherwise you get a site mismatch id
> - license checking requires a connection to the internet, which is not available on the compute nodes
> - number of cores allowed on the license is limited
> - number of simultaneous gurobi instances is limited

> **Snakemake**
> - shell scripts not running as expected (need to delete .snakemake everytime)