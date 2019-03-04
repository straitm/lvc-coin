#!/bin/bash

# Need to fail if nova fails in 'nova | tee'
set -o pipefail

# I know, I know.  Lots of positional arguments.  Not good.
infile=$1
fcl=$2
reco=$3
hist=$4
log=$5

if ! nova $infile -c $fcl -o $reco -T $hist 2> /dev/stdout | tee $log; then
  printf failed to process $infile
  exit 1
fi
