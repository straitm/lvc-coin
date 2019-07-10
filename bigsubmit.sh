#!/bin/bash

throttle()
{
  while [ $(ps f -u mstrait | wc -l) -gt 900 ]; do
    printf "Waiting for $(ps f -u mstrait | wc -l) to be less than 900\n"
    sleep 1m
  done
}

for e in GW150914 GW151012 GW151226 GW170104 GW170608 GW170729 GW170809 GW170814 GW170817 GW170818 GW170823 S190412m S190421ar S190425z S190426c S190503bf S190510g S190512at S190513bm S190517h S190519bj S190521g S190521r; do
  cat $SRT_PRIVATE_CONTEXT/ligo/goodnd.txt | while read m d; do
    throttle
    cd /nova/ana/users/mstrait/ligobgresults-$e
    if ! $SRT_PRIVATE_CONTEXT/ligo/should.sh $m $d neardet-ddactivity1; then
      continue
    fi
    getandredoloop.sh $m $d neardet-ddactivity1 2019 $e > /dev/null &
    echo Started $m $d near $e
    sleep 1
  done

  cat $SRT_PRIVATE_CONTEXT/ligo/goodfd.txt | while read m d; do
    cd /nova/ana/users/mstrait/ligobgresults-$e
    if ! $SRT_PRIVATE_CONTEXT/ligo/should.sh $m $d neardet-ddactivity1; then
      continue
    fi
    throttle
    getandredoloop.sh $m $d fardet-t02 2019 $e > /dev/null &
    echo Started $m $d far $e
    sleep 1
  done
done
