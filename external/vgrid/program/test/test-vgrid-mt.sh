#! /bin/sh 



n=2;
max=1;
min=1;
input="";
geom=two-halfs;
scale="1.0"

while [ $# -gt 0 ]
do
  case "$1" in
   -in) opt="$opt $1"
        input="$2"
        shift;;
   -out) opt="$opt $1"
         output="$2"
         shift;;
   -min) opt="$opt $1"
            min="$2"
            shift;;
   -max)  opt="$opt $1"
            max="$2"
            shift;;
   -n)    opt="$opt $1"
          n="$2"
          shift;;
   -geom) opt="$opt $1"
          geom="$2"
          shift;;
    -scale) opt="$opt $1"
           scale="$2"
           shift;;
     *) echo "$0: Unknown option $1"
        exit 1;;
  esac
  shift
done

if test "$input" = ""
then
  input="${n}x${n}x${n}-min${min}"
fi

mkdir -p out-mt
gridbase=`basename ${input} .v`
img_b=${gridbase}-${geom}.v
img=out-mt/${img_b}
grid_b=${img_b}-mesh-mt
grid=out-mt/${grid_b}
vgridout=out-mt/${grid_b}-olevel3.out

echo "./genimg -n $n -out ${img} -geom ${geom} -scale $scale" && \
      ./genimg -n $n -out ${img} -geom ${geom} -scale $scale;
stt=$?

echo "../vgrid -in ${img} -out ${grid}.ascii -format ascii -min $min -max $max -olevel 3 -smooth marching > ${vgridout} 2>&1" \&& \
      ../vgrid -in ${img} -out ${grid}.ascii -format ascii -min $min -max $max -olevel 3 -smooth marching > ${vgridout} 2>&1;
stt=`expr $stt + $?`

echo "../vgrid -in ${img} -out ${grid}.gmv -format gmv -min $min -max $max -olevel 3 -smooth marching > /dev/null 2>&1" \&& \
      ../vgrid -in ${img} -out ${grid}.gmv -format gmv -min $min -max $max -olevel 3 -smooth marching > /dev/null 2>&1;
stt=`expr $stt + $?`


echo "diff ${vgridout} exp/${vgridout}" \
 &&   diff ${vgridout} exp/${vgridout};
stt=`expr $stt + $?`

echo "diff ${grid}.ascii   exp/${grid}.ascii" && \
      diff ${grid}.ascii   exp/${grid}.ascii;
stt=`expr $stt + $?`


exit $stt

