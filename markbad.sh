#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

name=$(printf "$1" | sed 's@/$@@')

if ! [ -d $outhistdir/$name ]; then
  echo No hist directory for $name, doing nothing
  exit 1
fi

touch $outhistdir/$name/bad

rm -fv $outhadddir/${name}*
