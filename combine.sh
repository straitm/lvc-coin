#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

PATH=$PATH:$SRT_PRIVATE_CONTEXT/ligo

pnfsdir="$1"
basepnfsdir="$(basename "$1")"

if ! [ -e $pnfsdir/complete ] && ! iscomplete.sh $pnfsdir; then
  echo $pnfsdir is not complete
  exit 1
fi

if hadd  -f $outhadddir/$basepnfsdir.hadded.root $pnfsdir/*.root &> /dev/null; then
  echo $pnfsdir concatenated
else
  echo $pnfsdir concatenatation failed
  exit 1
fi
