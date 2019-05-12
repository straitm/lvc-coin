#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

PATH=$PATH:$SRT_PRIVATE_CONTEXT/ligo

pnfsdir="$1"
basepnfsdir="$(basename "$1")"

# Repair the timestamp in the name of the directory by replacing the two
# hyphens in the time-of-day with colons.  (They can't be colons in pnfs
# because of a SAM bug.)
basepnfsdir=${basepnfsdir/-/:}
basepnfsdir=${basepnfsdir/-/:}
basepnfsdir=${basepnfsdir/-/:}
basepnfsdir=${basepnfsdir/-/:}
basepnfsdir=${basepnfsdir/:/-}
basepnfsdir=${basepnfsdir/:/-}


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
