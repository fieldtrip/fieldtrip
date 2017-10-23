#!/bin/bash
set -e -u

# for cross-compilation you can use something like this
# MAKE="make $1 MACHINE=i386 PLATFORM=Darwin"

MAKE="make $1"
PLATFORM=`gcc -dumpmachine`
UNAME=`uname`
MACHINE=`uname -m`

echo $UNAME
echo $MACHINE

function contains () {
    # helper function to determine whether a bash array contains a certain element
    # http://stackoverflow.com/questions/3685970/bash-check-if-an-array-contains-a-value
    local n=$#
    local value=${!n}
    for ((i=1;i < $#;i++)) {
        if [ "${!i}" == "${value}" ]; then
            echo "y"
            return 0
        fi
    }
    echo "n"
    return 1
}

if [ "$UNAME" = "Linux" ]; then
  if [ "$MACHINE" = "armv6l" ]; then
    BLACKLIST=(amp audio biosemi ctf emotiv gtec neuralynx neuromag siemens tmsi tobi)
  elif [ "$MACHINE" = "armv7l" ]; then
    BLACKLIST=(amp audio biosemi ctf emotiv gtec neuralynx neuromag siemens tmsi tobi)
  elif [ "$MACHINE" = "x86_64" ]; then
    BLACKLIST=(audio ctf emotiv neuralynx siemens tmsi tobi)
  else
    BLACKLIST=(audio emotiv neuralynx siemens tmsi tobi)
  fi
fi

if [ "$UNAME" = "Darwin" ]; then
  BLACKLIST=(emotiv neuralynx siemens neuromag tmsi tobi ctf)
fi

if [ "$UNAME" = "MINGW32_NT-6.1" ]; then
  BLACKLIST=(amp audio emotiv siemens neuralynx neuromag tmsi tobi ctf)
fi

if [ "$UNAME" = "MINGW64_NT-6.1" ]; then
  BLACKLIST=(amp audio emotiv siemens neuralynx neuromag tmsi tobi ctf)
fi

echo Building buffer and ODM...
(cd buffer/src && $MAKE)
(cd buffer/cpp && $MAKE)

echo Building acquisition...
for ac in `ls -d acquisition/*/`; do
  SHORTNAME=`basename $ac`
  echo $SHORTNAME
  if [ $(contains ${BLACKLIST[@]} ${SHORTNAME}) == "y" ] ; then
    echo \'$SHORTNAME\' is blacklisted for this platform.
  else
    echo Building \'$ac\'...
    (cd $ac && $MAKE)
  fi
done

echo Building utilities...
for ac in `ls -d utilities/*/`; do
  SHORTNAME=`basename $ac`
  echo $SHORTNAME
  if [ $(contains ${BLACKLIST[@]} ${SHORTNAME}) == "y" ] ; then
    echo \'$SHORTNAME\' is blacklisted for this platform.
  else
    echo Building \'$ac\'...
    (cd $ac && $MAKE)
  fi
done
