#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

name=$(printf "$1" | sed 's@/$@@')

touch $outhistdir/$name/bad

rm -fv $outhadddir/${name}*
