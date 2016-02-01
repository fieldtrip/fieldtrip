#!/bin/bash
set -e -u

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
    BLACKLIST=(amp audio biosemi ctf emotiv neuralynx neuromag siemens tmsi tobi)
  else
    BLACKLIST=(audio emotiv neuralynx siemens tmsi tobi)
  fi
fi

if [ "$UNAME" = "Darwin" ]; then
  BLACKLIST=(audio emotiv neuralynx siemens neuromag tmsi tobi ctf)
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

