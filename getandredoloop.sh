#!/bin/bash

# Cache files, then run jobs to process them all, running redo loops that
# look for failed jobs and resubmit them. Try to do the redo loops gently
# enough that we don't fill up /tmp with overlapping submissions too often or
# cause other problems.

cd /nova/ana/users/mstrait/ligometalog/

if [ $# -eq 5 ] || [ $# -eq 6 ]; then
  month=$1
  day=$2
  trigger=$3
  year=$4
  export GWNAME=$5
  time=$6
  if [ $time ]; then
    unixtime=$(TZ=UTC date +%s -d"$month $day $time $year")
    fracsec=$(cut -d. -f 2 -s <<< $time)
    unixtime+=.$fracsec
  else
    unixtime=$(date +%s -d"$month $day 8:29:01 $year")
  fi
  rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$unixtime).${fracsec}Z
  . $SRT_PRIVATE_CONTEXT/ligo/env.sh
  log=$GWNAME-$month-$day-$year-$unixtime-$trigger.$(date "+%Y-%m-%dT%H:%M:%S").log
elif [ $# -eq 2 ]; then
  trigger=$1
  export GWNAME=$2
  export REALGWEVENT=1
  . $SRT_PRIVATE_CONTEXT/ligo/env.sh
  rfctime=$realgweventtime
  fracsec=$(cut -d. -f 2 -s <<< $realgweventtime | tr -d Z)
  unixtime=$(TZ=UTC date +%s -d "${realgweventtime/T/ }").$fracsec
  log=real-$GWNAME-$trigger.$(date "+%Y-%m-%dT%H:%M:%S").log
else
  echo Syntax: $(basename $0) month day trigger year gwname '[UTC time]'
  echo "For example: $(basename $0) Jan 1 fardet-ddsnews 2019 S190426c"
  echo "  gwname defines the sky map and effective time for pointing"
  echo "If UTC time is not given, the 8:30 SNEWS trigger time is used"
  echo "That's for background time samples"
  echo
  echo "Or syntax for real events: $(basename $0) trigger gwname"
  exit 1
fi

echo $(basename $0) running on $HOSTNAME

(
$SRT_PRIVATE_CONTEXT/ligo/getfilelist.sh $unixtime $trigger
ret=$?
if [ $ret -eq 1 ]; then
  echo getfilelist $unixtime $trigger failed
  exit 1
elif [ $ret -eq 2 ]; then
  # if getfilelist returned 2, it means it found no files to process
  exit 0
fi

defbase=strait-ligo-coincidence
recodef=$defbase-reco-$unixtime-$trigger
rawdef=$defbase-artdaq-$unixtime-$trigger

# getfilelist.sh only allows good definitions to exist
if samweb list-definitions | grep -q $recodef; then
  def=$recodef
else
  def=$rawdef
fi

if $SRT_PRIVATE_CONTEXT/ligo/redoloop_ligo.sh $def; then
  rfctimesafeforsam=${rfctime//:/-}
  if [ $trigger != neardet-t00 ]; then
    $SRT_PRIVATE_CONTEXT/ligo/combine.sh \
      $outhistdir/$rfctimesafeforsam-$trigger
  fi
fi
) 2> /dev/stdout | tee /nova/ana/users/mstrait/ligometalog/$log
