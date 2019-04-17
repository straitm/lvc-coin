#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

PATH=$PATH:$SRT_PRIVATE_CONTEXT/ligo

pnfsdir="$1"
basepnfsdir="$(basename "$1")"

if ! [ -e $pnfsdir/complete ] && ! iscomplete.sh $pnfsdir; then
  echo $pnfsdir is not complete
  exit 1
fi

hadd -f $outhadddir/$basepnfsdir.hadded.root $pnfsdir/*.root
