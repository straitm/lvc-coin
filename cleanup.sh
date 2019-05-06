#!/bin/bash

# Removes duplicate output files in the current directory that differ only
# by the release name.

for thing in hists reco; do

  if ! ls *det*.$thing.root &> /dev/null; then
    exit
  fi

  ls *det*.$thing.root | while read f; do
    echo $f | cut -d_ -f1-4
  done | sort | uniq -c | while read n base; do
    ls ${base}*.$thing.root | tail -n +2
  done | xargs rm -fv

done
