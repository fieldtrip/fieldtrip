#!/bin/sh
#
# shell script to start a peer in slave mode
# it immediately restart if Matlab exits

while [ true ] ; do 
  /opt/matlab78/bin/matlab -nodesktop -r "cd /home/common/matlab/fieldtrip/peer; peerslave('threads', 1);"
done

