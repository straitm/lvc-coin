#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

PATH=$PATH:$SRT_PRIVATE_CONTEXT/ligo

pnfsdir="$1"
basepnfsdir="$(basename "$1")"

if ! [ -e $pnfsdir/complete ] && ! iscomplete.sh $pnfsdir; then
  echo $pnfsdir is not complete
  exit 1
fi

if [ -e $pnfsdir/bad ]; then
  echo $pnfsdir is marked bad, skipping
  exit 0
fi

out=$outhadddir/$basepnfsdir.hadded.root 

mkdir -p $outhadddir

if [ -e $out ]; then
  echo $out already exists, skipping
  exit 0
fi

(cd $pnfsdir; cleanup.sh)

if timeout 120 hadd  -f $out $pnfsdir/*hists.root; then
  echo $pnfsdir concatenated
else
  echo $pnfsdir concatenatation failed
  rm -fv $out
  exit 1
fi
