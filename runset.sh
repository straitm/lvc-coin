#!/bin/bash

if ! [ $1 ]; then
  echo specify an analysis type
  echo \(That is, the part of the name of the fcl between ligojob_ and .fcl\)
  exit 1
fi

type=$1

shift

for f in $@; do
  base=$(basename $f .artdaq.root)
  nova -c ligojob_$type.fcl $f \
       -o $base.ligo.root \
       -T $base.hists.root
done
