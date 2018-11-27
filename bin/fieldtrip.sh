#!/bin/sh
#
# Script for execution of deployed FieldTrip applications. It sets
# up the MCR environment for the current $ARCH and executes the
# specified command.
#
# Use as
#   fieldtrip.sh <MATLABROOT> script.m
#   fieldtrip.sh <MATLABROOT> script1.m script2.m ...
#   fieldtrip.sh <MATLABROOT> jobfile.mat
#
# See also the ft_standalone MATLAB function

# Copyright (C) 2011-2018, Robert Oostenveld
#
# This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
# for the documentation and details.
#
#    FieldTrip is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    FieldTrip is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
#

exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  echo Setting up environment variables
  MCRROOT="$1"
  echo ---
  DYLD_LIBRARY_PATH=.:${MCRROOT}/runtime/maci64 ;
  DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${MCRROOT}/bin/maci64 ;
  DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${MCRROOT}/sys/os/maci64;
  XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
  export DYLD_LIBRARY_PATH;
  export XAPPLRESDIR;
  echo DYLD_LIBRARY_PATH is ${DYLD_LIBRARY_PATH};
  shift 1
  args=
  while [ $# -gt 0 ]; do
      token=$1
      args="${args} ${token}" 
      shift
  done
  "${exe_dir}"/fieldtrip.app/Contents/MacOS/fieldtrip $args
fi
exit

