#!/bin/bash

if ! [ $1 ]; then
  echo "$(basename $0) [type] [time] [files]"
  echo "  type:  The analysis type, that is, the part of"
  echo "         the name of the fcl between ligojob_ and .fcl"
  echo "  time:  Unix time stamp of the GW event or whatever"
  echo "  files: Any number of art files"
  exit 1
fi

# Should only be needed temporarily
CET_PLUGIN_PATH+=:$SRT_PRIVATE_CONTEXT/lib/$SRT_SUBDIR
export CET_PLUGIN_PATH

type=$1
shift

if ! [ $1 ]; then
  echo specify the unix time of the event
  exit 1
fi
unixtime=$1
shift

# CHANGE THIS to the appropriate skymap
# TODO: should be a command line argument
skymap=/nova/ana/users/mstrait/skymaps/LALInference_skymap-GW170817.fits

export TZ=UTC
fracsec=$(cut -d. -f 2 <<< $unixtime)
rfctime=$(date "+%Y-%m-%dT%H:%M:%S" -d @$unixtime).${fracsec}Z

if [ $? -ne 0 ]; then
  exit 1
fi

if ! $SRT_PRIVATE_CONTEXT/ligo/makefcl.sh $type $unixtime $skymap; then
  exit 1
fi

fcl=ligojob_$type.$rfctime.fcl

if ! [ $1 ]; then
  echo You did not specify any files to processes.  Technical success\!
  exit 0
fi

out=$type.ligohists.root

for f in $@; do
  if ! [ -e "$f" ]; then
    echo $f does not exist > /dev/stderr
    exit 1
  fi
done

for f in $@; do
  base=$(basename $(basename $(basename $f .root) .raw) .artdaq)
  reco=$base.ligo${type}.root
  if echo $type | grep -q noreco; then
    reco=/dev/null
  fi
  log=$base.ligo.$type.log
  hist=$base.hists.root
  hists="$hists $hist"

  if [ -e $hist ]; then
    echo $hist already exists. Skipping $base. Will still include in sum. > /dev/stderr
    continue
  fi

  if [ -e $log ]; then
    echo $log already exists. Skipping $base. Will still include in sum. > /dev/stderr
    continue
  fi

  echo $SRT_PRIVATE_CONTEXT/ligo/runone.sh $f $fcl $reco $hist $log
done | ~mstrait/bin/parallel -j 3 --line-buffer

if [ $? -eq 0 ]; then
  hadd -f $out $hists
  echo Full histograms in $out
fi
