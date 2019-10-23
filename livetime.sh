#!/bin/bash

(root -b $1 << EOF
  ligoanalysis->cd()
  blindlive->Integral()
  .q
EOF
) 2> /dev/null | tail -n 1 | awk '{print $2}'
