#!/bin/sh
#
# helper script to start multiple peers in slave mode

command=$1

# 512MB =  536870912
# 1GB   = 1073741824
# 2GB   = 2147483648 
# 4GB   = 4294967296 
# 8GB   = 8589934592 

case "$HOSTNAME" in 
mentat001|mentat002)
    # has 16GB of RAM and 4 cores
	peermem[0]="auto"
	peermem[1]="auto"
	peermem[2]="auto"
	peermem[3]="auto"
	;;
mentat003|mentat004)
    # has 32GB of RAM and 8 cores
	peermem[0]="auto"
	peermem[1]="auto"
	peermem[2]="auto"
	peermem[3]="auto"
	peermem[4]="auto"
	peermem[5]="auto"
	peermem[6]="auto"
	peermem[7]="auto"
	;;
mentat005)
    # has 48GB of RAM and 12 cores
	peermem[ 0]=8589934592
	peermem[ 1]=8589934592
	peermem[ 2]=4294967296
	peermem[ 3]=4294967296
	peermem[ 4]=4294967296
	peermem[ 5]=4294967296
	peermem[ 6]="auto"
	peermem[ 7]="auto"
	peermem[ 8]="auto"
	peermem[ 9]="auto"
	peermem[10]="auto"
	peermem[11]="auto"
	;;
esac

##############################################################################

PLATFORM=`gcc -dumpmachine`
HOSTNAME=`hostname -s`

case "$PLATFORM" in
"x86_64-redhat-linux")
    MATLABARCH=glnxa64
	MATLABPATH=/opt/matlab77
	;;
"i386-redhad-linux")
    MATLABARCH=glnx86
	MATLABPATH=/opt/matlab77
	;;
"i686-apple-darwin9")
    MATLABARCH=maci
	MATLABPATH=/opt/matlab77
	;;
"i686-apple-darwin10")
    MATLABARCH=maci
	MATLABPATH=/opt/matlab77
	;;
*)
    MATLABARCH=unknown
	MATLABPATH=unknown
	;;
esac

export LD_LIBRARY_PATH=$MATLABPATH/bin/$MATLABARCH/ 
EXECPATH=/home/common/matlab/fieldtrip/peer
EXECUTABLE=$EXECPATH/peerslave.$MATLABARCH
OPTIONS=

##############################################################################

case $command in
start)
length=${#peermem[@]}
index=0
while [ "$index" -lt "$length" ] ; do
	memavail=${peermem[$index]}
	if [[ "$memavail" -eq "auto" ]] ; then
		OPTIONS="--smartmem 1"
	else
		OPTIONS="--smartmem 0 --memavail $memavail"
	fi
	eval $EXECUTABLE $OPTIONS &
	peerstatus[$index]=$?
	peerpid[$index]=$!
	let "index = $index + 1"
done
#echo ${peermem[@]}
#echo ${peerpid[@]}
#echo ${peerstatus[@]}
;;

stop)
killall $EXECUTABLE
;;

*)
echo Use as $0 "<start|stop>"
;;
esac

