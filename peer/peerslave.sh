#!/bin/sh
#
# shell script to start a peer in slave mode, restart if it exits

while [ true ] ; do 
  /opt/matlab78/bin/matlab -nodesktop -r "peerslave('threads', 1)"
done

