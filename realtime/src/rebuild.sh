#!/bin/bash
set -e -u
MAKE="make $1"
PLATFORM=`uname`

echo Building buffer and ODM...
(cd buffer/src && $MAKE)
(cd buffer/cpp && $MAKE)

echo Building acquisition software...
if [[ $PLATFORM == 'Linux' ]]; then
  BLACKLIST='audio emotiv neuralynx siemens tmsi tobi'
fi

if [[ $PLATFORM == 'Darwin' ]]; then
  BLACKLIST='audio emotiv neuralynx siemens neuromag tmsi tobi'
fi

#for ac in `ls acquisition`; do
for ac in `ls -d acquisition/*/`; do
  if [[ $BLACKLIST ==  *`basename $ac`* ]]; then
    echo \'$ac\' is blacklisted for this platform.
  else
    echo Building \'$ac\'...
    (cd $ac && $MAKE)
  fi
done;

# TODO: add utilities
