#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

markcomplete.sh

cd $outhistdir

for f in */complete; do
  $SRT_PRIVATE_CONTEXT/ligo/combine.sh $(dirname $f)
done
