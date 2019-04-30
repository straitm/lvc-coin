#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

if [ $# -ne 1 ] && [ $# -ne 2 ]; then
  echo Syntax: $(basename $0) unixtime {def}
  echo Override the definition derived from unixtime with def.
  exit 1
fi

unixtime=$1

if [ $2 ]; then
  def=$2
fi

def=strait-ligo-coincidence-artdaq-${unixtime}-neardet-t00

rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$unixtime).${fracsec}Z
fcl=spilltimejob.fcl

filesperjob=$(samweb list-files defname: $def | wc -l)
lifetime=$(( 300 * $filesperjob ))

tag=$(cat $SRT_PRIVATE_CONTEXT/.base_release)
testrel=/nova/app/users/mstrait/novasoft-ligo/

# This is needed because the way submission works these days is that it
# tars up the whole test release in /tmp, which only has 2GB.
fails=0
while true; do
  kfree=$(df /tmp | tail -n 1 | awk '{print $3}')
  if [ $kfree -lt 1000000 ]; then
    echo /tmp is too full.  Waiting for it to clear out
    df -h /tmp | tail -n 1
    sleep 1m
  elif ls -l /tmp | grep 'mstrait.*tmp......$'; then
    echo waiting for another submit to finish dumping in /tmp
    sleep 1m
  else
    break
  fi
  let fails++
  if [ $fails -gt 100 ]; then
    echo OK, I have had it.  Going to clear my files out of /tmp
    rm -fv /tmp/tmp* /tmp/mstrait.* /tmp/joblist.* /tmp/headtail.*
    fails=0
  fi
done

# Remember to take out "--test_submission" and/or --nevts and/or set --njobs to
# $njobs

submit_nova_art.py \
--opportunistic \
--maxopt \
--logs \
--copyOutScript /bin/true \
--jobname $def \
--defname $def \
--njobs 1 \
--files_per_job $filesperjob \
--tag $tag \
--testrel $testrel \
--config $fcl \
--expected_lifetime $lifetime
