#!/bin/bash

if ! [ -e "$1" ]; then
  echo "$1": File not found
  exit 1
fi

(root -b $1 << EOF
  snanalysis->cd()
  supernovalikelive->Integral(-1, 1001)
  .q
EOF
) 2> /dev/null | tail -n 1 | awk '{print $2}'
