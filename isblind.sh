#!/bin/bash

if ! [ $1 ]; then
  echo give me a file name
  exit 1
fi

printf 'ligoanalysis->ls()\n.q\n' | root -b -n -l "$1" | grep -q blind
