#!/bin/bash

if [ $# -ne 1 ] && [ $# -ne 2 ]; then
  echo Syntax: $(basename $0) sam_definition {skymap}
  echo I will use a dummy skymap if you dont give one
  exit 1
fi

def=$1
skymap=$2

unixtime=$(printf $def | cut -d- -f 5)
analysis_type_key=$(printf $def | cut -d- -f 6-7)

if [ $analysis_type_key == fardet-t02 ]; then
  type=minbiasfd
elif [ $analysis_type_key == fardet-ddsnews ]; then
  type=minbiasfd
elif [ $analysis_type_key == fardet-ddenergy ]; then
  type=ddenergy
elif [ $analysis_type_key == neardet-ddactivity1 ]; then
  type=ndactivity
elif [ $analysis_type_key == neardet-ddsnews ]; then
  type=minbiasnd
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

# Remember to take out "--test_submission" and/or --nevts and/or set --njobs to
# $njobs

submit_nova_art.py \
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
--dest /pnfs/nova/scratch/users/mstrait/ligo/ \
--expected_lifetime 7200
