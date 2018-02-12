#!/bin/bash

if ! [ $1 ]; then
  echo specify an analysis type
  echo \(That is, the part of the name of the fcl between ligojob_ and .fcl\)
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

fcl=ligojob_$type.$rfctime.fcl

if ! [ -e $SRT_PRIVATE_CONTEXT/job/ligojob_$type.fcl ]; then
  echo Analysis type $type not supported
  exit 1
fi

cat $SRT_PRIVATE_CONTEXT/job/ligojob_$type.fcl | \
  sed "/this_here_ligo: @local/a this_here_ligo.GWEventTime: \"$rfctime\"" > $fcl

# Need to fail if nova fails in 'nova | tee'
set -o pipefail

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
  if ! nova $f -c $fcl $recoopt -T $hist 2> /dev/stdout | tee $log; then
    exit 1
  fi
done

out=$type.ligohists.root

hadd -f $out $hists
echo Full histograms in $out
