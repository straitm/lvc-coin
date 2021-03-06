#!/bin/bash

# Cache files, then run jobs to process them all, running redo loops that
# look for failed jobs and resubmit them. Try to do the redo loops gently
# enough that we don't fill up /tmp with overlapping submissions too often or
# cause other problems.

cd /nova/ana/users/mstrait/ligometalog/

hourbefore()
{
  fracsec=$(echo $1 | cut -d. -f 2 | cut -dZ -f 1)
  indate=$(echo $1 | cut -d. -f 1 | sed 's/T/ /')
  export TZ=UTC
  printf "%s.%sZ\n" \
    $(date "+%Y-%m-%dT%H:%M:%S" -d @$(($(date +%s -d "$indate UTC") - 3600))) \
    $fracsec
}


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
  logbase=$GWNAME-$month-$day-$year-$unixtime-$trigger
  log=$logbase.$(date "+%Y-%m-%dT%H:%M:%S").log
elif [ $# -eq 2 ]; then
  trigger=$1
  export GWNAME=$2
  export REALGWEVENT=1
  . $SRT_PRIVATE_CONTEXT/ligo/env.sh
  rfctime=$realgweventtime
  fracsec=$(cut -d. -f 2 -s <<< $realgweventtime | tr -d Z)
  unixtime=$(TZ=UTC date +%s -d "${realgweventtime/T/ }").$fracsec
  logbase=real-$GWNAME-$trigger
  log=$logbase.$(date "+%Y-%m-%dT%H:%M:%S").log
elif [ $# -eq 3 ]; then
  if [ $3 != side ]; then
    echo I expected the third argument to be literally \"side\"
    exit 1
  fi
  trigger=$1
  export GWNAME=$2
  export SIDEBAND=1
  . $SRT_PRIVATE_CONTEXT/ligo/env.sh
  rfctime=$(hourbefore $realgweventtime)
  fracsec=$(cut -d. -f 2 -s <<< $realgweventtime | tr -d Z)
  unixtime=$(TZ=UTC date +%s -d "${rfctime/T/ }").$fracsec
  logbase=sideband-$GWNAME-$trigger
  log=$logbase.$(date "+%Y-%m-%dT%H:%M:%S").log
else
  echo Syntax: $(basename $0) month day trigger year gwname '[UTC time]'
  echo "For example: $(basename $0) Jan 1 fardet-ddsnews 2019 S190426c"
  echo "  gwname defines the sky map and effective time for pointing"
  echo "If UTC time is not given, the 8:30 SNEWS trigger time is used"
  echo "That's for background time samples"
  echo
  echo "Or for real events: $(basename $0) trigger gwname"
  echo "For blind livetime: GWBLIND=1 $(basename $0) trigger gwname"
  echo "Or for sideband: $(basename $0) trigger gwname \"side\""
  exit 1
fi

echo $(basename $0) running on $HOSTNAME

if ! klist &> /dev/null; then
  echo You should renew your kerberos ticket
  exit 1
fi

if ps f -u mstrait | grep tee | grep -B 3 ligometalog/$logbase; then
  echo It looks like you are already running $(basename $0) with these
  echo options.  Trying to run multiple copies is, at best, wasteful, and
  echo might also be actively counterproductive, so I am quitting.
  exit 1
fi

if [ $trigger == neardet-t00 ] && [ -e $spilldir/spills-$unixtime-$rfctime.txt ]; then
  echo Looks like a spill file has already been generated for this time
  exit 0
fi

fulllog=/nova/ana/users/mstrait/ligometalog/$log
echo Output going to log file. You might want to do:
echo tail -f $fulllog

exec > $fulllog 2>&1

$SRT_PRIVATE_CONTEXT/ligo/getfilelist.sh $unixtime $trigger

ret=$?

if [ $ret -eq 2 ]; then
  # if getfilelist returned 2, it means it found no files to process
  exit 0
elif [ $ret -gt 0 ]; then
  echo getfilelist $unixtime $trigger failed
  exit 1
fi

defbase=strait-ligo-coincidence
recodef=$defbase-reco-$unixtime-$trigger-$gwbase
rawdef=$defbase-artdaq-$unixtime-$trigger

# getfilelist.sh only allows good definitions to exist
if samweb describe-definition $recodef &> /dev/null; then
  def=$recodef
else
  def=$rawdef
fi

if $SRT_PRIVATE_CONTEXT/ligo/redoloop_ligo.sh $def $rawdef; then
  rfctimesafeforsam=${rfctime//:/-}
  if [ $trigger != neardet-t00 ]; then
    $SRT_PRIVATE_CONTEXT/ligo/combine.sh \
      $outhistdir/$rfctimesafeforsam-$trigger
  fi
fi
