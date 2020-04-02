#!/bin/bash

# override system install
export PATH=/usr/local/bin:$PATH

# install deps. note than many are included as a part of brew-core
# these days. hence this list isn't comprehensive
brew update
brew install mpich swig svn udunits
# matplotlib currently doesn't have a formula
# teca fails to locate mpi4py installed from brew
pip3 install numpy mpi4py matplotlib torch

# install data files.
# On Apple svn is very very slow. On my mac book pro
# this command takes over 14 minutes, while on a Linux
# system on the same network it takes less than 2 minutes.
# travis will kill a build  if it does not get console output
# for 10 min. The following snippet sends progress marks
# to the console while svn runs.
echo 'svn co svn://svn.code.sf.net/p/teca/TECA_data@${TECA_DATA_REVISION} TECA_data &'
svn svn co svn://svn.code.sf.net/p/teca/TECA_data@${TECA_DATA_REVISION} TECA_data &
svn_pid=$!
while [ -n "$(ps -p $svn_pid -o pid=)" ]
do
  echo -n "."
  sleep 2s
done
