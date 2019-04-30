#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

name=$(printf "$1" | sed 's@/$@@')
detdate=$(printf "$name" | cut -d- -f 1-4)

if echo $name | grep -q neardet; then
  for stream in ddactivity1 ddsnews; do 
    d=$outhistdir/${detdate}-$stream
    set -x
    mkdir -p $d
    touch $d/bad
    set +x
  done
else
  for stream in t02 ddenergy ddsnews; do 
    d=$outhistdir/${detdate}-$stream
    set -x
    mkdir -p $d
    touch $d/bad
    set +x
  done
fi

rm -fv $outhadddir/${detdate}*
