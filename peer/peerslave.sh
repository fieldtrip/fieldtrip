#!/bin/sh
#
# shell script to start a peer in slave mode
# it restarts if a problem is detected

while [ true ] ; do 
  /opt/matlab78/bin/matlab -nodesktop -r "addpath('/home/common/matlab/fieldtrip/peer'); try, peerslave('threads', 1, 'maxidle', 600); end; "
done

