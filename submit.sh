#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

if [ $# -ne 2 ] && [ $# -ne 3 ] && [ $# -ne 4 ]; then
  echo Syntax: $(basename $0) unixtime analysis_type {skymap} {def}
  echo Where analysis_type is fardet-t02, fardet-ddsnews, fardet-ddenergy
  echo neardet-ddactivity1 or neardet-ddsnews.
  echo I will use a dummy skymap if you do not give one
  echo Override the definition derived from unixtime with def.
  exit 1
fi

unixtime=$1
analysis_type_key=$2

def=strait-ligo-coincidence-artdaq-${unixtime}-${analysis_type_key}
skymap=/pnfs/nova/users/mstrait/ligo/LALInference_skymap-GW170817.fits
if [ $3 ]; then
  skymap=$3
fi
if [ $4 ]; then
  def=$4
fi

if [ $analysis_type_key == fardet-t02 ]; then
  type=minbiasfd
  lifetime=40000
elif [ $analysis_type_key == fardet-ddsnews ]; then
  type=minbiasfd_rawinput
  lifetime=172800
elif [ $analysis_type_key == fardet-ddenergy ]; then
  type=ddenergy
  lifetime=14400
elif [ $analysis_type_key == neardet-ddactivity1 ]; then
  type=ndactivity
  lifetime=14400
elif [ $analysis_type_key == neardet-ddsnews ]; then
  type=minbiasnd
  lifetime=14400
else
  echo I cannot figure out what analysis type to run for $analysis_type_key
  echo which I got from $def
  exit 1
fi

# Tell the fcl we are reading the skymap from the CWD since it will be shipped
# there by the grid.  (*Don't* do that for interactive jobs.)
if ! $SRT_PRIVATE_CONTEXT/ligo/makefcl.sh $type $unixtime $(basename $skymap); then
  exit 1
fi

rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$unixtime).${fracsec}Z
fcl=ligojob_$type.$rfctime.fcl

njobs=$(samweb list-files defname: $def | wc -l)
tag=$(cat $SRT_PRIVATE_CONTEXT/.base_release)
testrel=/nova/app/users/mstrait/novasoft-ligo/

outdir=$outhistdir/$rfctime-$analysis_type_key

mkdir -p $outdir

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
  fi
done

# Remember to take out "--test_submission" and/or --nevts and/or set --njobs to
# $njobs

submit_nova_art.py \
--inputfile $skymap \
--opportunistic \
--maxopt \
--logs \
--histTier hists \
--copyOut \
--jobname $def \
--defname $def \
--njobs $njobs \
--files_per_job 1 \
--tag $tag \
--testrel $testrel \
--config $fcl \
--dest $outdir \
--expected_lifetime $lifetime
