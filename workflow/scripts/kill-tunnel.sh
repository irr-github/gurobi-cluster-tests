PORT=$1
SSH_TUNNEL_PID=$(lsof -t -i :$PORT)l
kill $SSH_TUNNEL_PID
