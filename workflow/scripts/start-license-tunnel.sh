#set up the tunnel
PORT=$1
ssh -fN -D $PORT $USER@login01

# Find the PID of the tunnel using the specified port
SSH_TUNNEL_PID=$(lsof -t -i :$PORT)l
export https_proxy=socks5://127.0.0.1:$PORT
export SSL_CERT_FILE=/p/projects/rd3mod/ssl/ca-bundle.pem_2022-02-08
export GRB_CAFILE=/p/projects/rd3mod/ssl/ca-bundle.pem_2022-02-08

# start gurobi script as before
export GUROBI_HOME="/p/projects/rd3mod/gurobi1103/linux64"
export PATH="$PATH:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
export GRB_LICENSE_FILE=/p/projects/rd3mod/gurobi_rc/gurobi.lic
#export GRB_LICENSE_FILE=/home/ivanra/gurobi.lic
export GRB_CURLVERBOSE=1

PATH=$GUROBI_HOME/bin:$PATH;export PATH
LD_LIBRARY_PATH=$GUROBI_HOME/lib:$LD_LIBRARY_PATH;export LD_LIBRARY_PATH
PYTHONHOME=$GUROBI_HOME;export PYTHONHOME
PYTHONPATH=$GUROBI_HOME:$PYTHONPATH;export PYTHONPATH

PYTHONSTARTUP=$PYTHONHOME/lib/gurobi.py;export PYTHONSTARTUP

#$PYTHONHOME/bin/python3.7 "$@"
# Use the currently active Python interpreter (from `which python`)
$(which python) "$@"
echo "setup"
