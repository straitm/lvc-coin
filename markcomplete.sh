#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

cd $outhistdir

for f in *; do
  if ! [ -d $f ]; then
    continue
  fi
  if ! [ -e $f/bad ] && ! [ -e $f/complete ] && $SRT_PRIVATE_CONTEXT/ligo/iscomplete.sh $PWD/$f; then
    touch $f/complete
  fi
done
