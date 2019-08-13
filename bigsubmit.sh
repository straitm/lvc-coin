#!/bin/bash

throttle() {
  proc=$(ps f -u mstrait | wc -l)
  load=$(cat /proc/loadavg | cut -d. -f1)
  while [ $proc -ge 700 ] || [ $load -gt 2 ]; do
    printf "Waiting for $proc to be less than 700 and $load less than 3\n"
    sleep 1m
    proc=$(ps f -u mstrait | wc -l)
    load=$(cat /proc/loadavg | cut -d. -f1)
  done
}


for e in S190602aq S190630ag S190701ah S190706ai S190707q; do
  cat $SRT_PRIVATE_CONTEXT/ligo/goodnd.txt | while read m d; do
    throttle
    cd /nova/ana/users/mstrait/ligobgresults-$e
    getandredoloop.sh $m $d neardet-ddactivity1 2019 $e > /dev/null &
    echo Started $m $d near $e
    sleep 1
  done

  cat $SRT_PRIVATE_CONTEXT/ligo/goodfd.txt | while read m d; do
    cd /nova/ana/users/mstrait/ligobgresults-$e
    throttle
    getandredoloop.sh $m $d fardet-t02 2019 $e > /dev/null &
    echo Started $m $d far $e
    sleep 1
  done
done
