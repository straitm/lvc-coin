#!/bin/bash

(root -b $1 << EOF
  snanalysis->cd()
  supernovalikelive->Integral(-1, 1001)
  .q
EOF
) 2> /dev/null | tail -n 1 | awk '{print $2}'
