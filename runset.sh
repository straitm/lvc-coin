#!/bin/bash

if ! [ $1 ]; then
  echo "$(basename $0) [type] [time] {novaargs} [files]"
  echo "  type:  The analysis type, that is, the part of"
  echo "         the name of the fcl between ligojob_ and .fcl"
  echo "  time:  Unix time stamp of the GW event or whatever"
  echo "  novaargs: optional extra arguments to the nova executable"
  echo "            must start with a \"-\""
  echo "  files: Any number of art files"
  echo "         must not start with a \"-\""
  exit 1
fi

type=$1
shift

if ! [ $1 ]; then
  echo specify the unix time of the event
  exit 1
fi
unixtime=$1
shift

export TZ=UTC
rfctime=$(date "+%Y-%m-%dT%H:%M:%SZ" -d @$unixtime)

# CHANGE THIS to the appropriate skymap
skymap=/nova/ana/users/mstrait/skymaps/LALInference_skymap-GW170817.fits

if [ $? -ne 0 ]; then
  exit 1
fi

fcl=ligojob_$type.$rfctime.fcl

if ! [ -e $SRT_PRIVATE_CONTEXT/job/ligojob_$type.fcl ]; then
  echo Analysis type $type not supported
  exit 1
fi

cat $SRT_PRIVATE_CONTEXT/job/ligojob_$type.fcl | \
  sed "/this_here_ligoanalysis: @local/a this_here_ligoanalysis.GWEventTime: \"$rfctime\"\
      \nthis_here_ligofilter.GWEventTime: \"$rfctime\"\
      \nthis_here_ligoanalysis.SkyMap: \"$skymap\"" > $fcl

# Need to fail if nova fails in 'nova | tee'
set -o pipefail

if [ "${1:0:1}" == - ]; then
  novaargs="$1"
  shift
fi

if ! [ $1 ]; then
  echo You did not specify any files to processes.  Success\!
  exit 0
fi

for f in $@; do
  if ! [ -e "$f" ]; then
    echo $f does not exist
    exit 1
  fi

  base=$(basename $f .artdaq.root)
  reco=$base.ligo.root
  if ! echo $type | grep -q noreco; then
    recoopt="-o $reco"
  fi
  log=$base.ligo.$type.log
  hist=$base.hists.root
  hists="$hists $hist"
  if ! nova $novaargs $f -c $fcl $recoopt -T $hist 2> /dev/stdout | tee $log; then
    exit 1
  fi
done

out=$type.ligohists.root

hadd -f $out $hists
echo Full histograms in $out
