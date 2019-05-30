#!/bin/bash

if [ $1 ]; then
  GWNAME=$1
elif [ $GWNAME ]; then
  :
else
  GWNAME=$(basename "$PWD")
  if [ $GWNAME == ${GWNAME##*-} ]; then
    echo Please set GWNAME
    exit 1
  fi
  export GWNAME=${GWNAME##*-}
  echo I got GWNAME = $GWNAME from the PWD.  I hope that is right.

  if basename "$PWD" | grep -q sideband; then
    export SIDEBAND=1
    echo Looks like this is a sideband sample.  I hope that is also right.
  elif basename "$PWD" | grep -qv bg; then
    export REALGWEVENT=1
    echo Looks like this is a real event.  I hope that is also right.
  fi
fi


. $SRT_PRIVATE_CONTEXT/ligo/env.sh

markcomplete.sh

mkdir -p $outhistdir
cd $outhistdir

for f in */complete; do
  $SRT_PRIVATE_CONTEXT/ligo/combine.sh $(dirname $f)
done
