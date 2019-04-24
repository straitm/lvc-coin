#!/bin/bash

# Removes duplicate output hist files in the current directory that differ only
# by the release name.

if ! ls *det*.root &> /dev/null; then
  exit
fi

ls *det*.root | while read f; do
  echo $f | cut -d_ -f1-4
done | sort | uniq -c | while read n base; do
  ls ${base}*root | tail -n +2
done | xargs rm -fv
