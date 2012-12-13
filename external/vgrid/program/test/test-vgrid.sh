#! /bin/sh 

mindflt=1

n=$1
min=${2} #-${mindflt}}
griddflt=${n}x${n}x${n}-min${min}
grid=${3} #-${griddflt}


if \
(echo "./genimg -n $n -out $grid.v" && ./genimg -n $n -out $grid.v \
&& \
echo "../vgrid -in ${grid}.v -out ${grid}mesh.ascii -format ascii -min $min -max $n -olevel 3 > ${grid}-olevel3.out" &&\
      ../vgrid -in ${grid}.v -out ${grid}mesh.ascii -format ascii -min $min -max $n -olevel 3 > ${grid}-olevel3.out 2>&1) \
then
echo "diff ${grid}-olevel3.out ${grid}-olevel3.out.exp" \
 && diff ${grid}-olevel3.out ${grid}-olevel3.out.exp;
stt=$?
echo "diff ${grid}mesh.ascii   ${grid}mesh.ascii.exp" && \
 diff ${grid}mesh.ascii   ${grid}mesh.ascii.exp;
stt=`expr $stt + $?`
exit $stt
else
 exit 1
fi
