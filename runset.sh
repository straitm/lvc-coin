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
CET_PLUGIN_PATH+=:$SRT_PRIVATE_CONTEXT/lib/Linux2.6-GCC-maxopt
export CET_PLUGIN_PATH

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
# TODO: should be a command line argument
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

if ! [ $1 ]; then
  echo You did not specify any files to processes.  Technical success\!
  exit 0
fi

for f in $@; do
  if ! [ -e "$f" ]; then
    echo $f does not exist > /dev/stderr
    exit 1
  fi

  base=$(basename $f .artdaq.root)
  reco=$base.ligo.root
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

out=$type.ligohists.root

hadd -f $out $hists
echo Full histograms in $out
