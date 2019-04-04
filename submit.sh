#!/bin/bash

if [ $# -ne 2 ] && [ $# -ne 3 ] && [ $# -ne 4 ]; then
  echo Syntax: $(basename $0) unixtime analysis_type {skymap} {def}
  echo Where analysis_type is fardet-t02, fardet-ddsnews, fardet-ddenergy
  echo neardet-ddactivity1 or neardet-ddsnews.
  echo I will use a dummy skymap if you dont give one
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
  lifetime=20000
elif [ $analysis_type_key == fardet-ddsnews ]; then
  type=minbiasfd_rawinput
  lifetime=86400
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

if ! $SRT_PRIVATE_CONTEXT/ligo/makefcl.sh $type $unixtime $skymap; then
  exit 1
fi

rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$unixtime).${fracsec}Z
fcl=ligojob_$type.$rfctime.fcl

njobs=$(samweb list-files defname: $def | wc -l)
tag=$(cat $SRT_PRIVATE_CONTEXT/.base_release)
testrel=/nova/app/users/mstrait/novasoft-ligo/

outdir=/pnfs/nova/scratch/users/mstrait/ligo/$rfctime-$analysis_type_key

mkdir -p $outdir

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
