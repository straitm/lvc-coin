#!/bin/bash

for f in $@; do
  base=$(basename $f .artdaq.root)
  if ! [ -e $base.bpf.root ]; then
    nova -c ligojob.fcl $f \
         -o $base.bpf.root \
         -T $base.hists.root
  fi
done
