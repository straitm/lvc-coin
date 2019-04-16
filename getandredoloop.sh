#!/bin/bash

# Cache the files needed for the 8:30 DDsnews trigger of the given month and
# day in 2019 (or another year if specified), then run jobs to process them
# all, running redo loops that look for failed jobs and resubmit them.  Try to
# do the redo loops gently enough that we don't fill up /tmp with overlapping
# submissions too often or cause other problems.

cd /nova/ana/users/mstrait/ligometalog/

if [ $# -lt 3 ]; then
  echo Syntax: $(basename $0) month day trigger '[year]'
  echo "For example: $(basename $0) Jan 1 fardet-ddsnews"
  echo "Or:          $(basename $0) Oct 9 fardet-ddsnews 2018"
  exit 1
elif [ $# -gt 4 ]; then
  echo $(basename $0): Warning, ignoring extraneous arguments after the fourth
fi

month=$1
day=$2
trigger=$3
if [ $4 ]; then
  year=$4
else
  year=2019
fi

unixtime=$(date +%s -d"$month $day 8:29:01 $year")

(
if ! $SRT_PRIVATE_CONTEXT/ligo/getfilelist.sh $unixtime $trigger; then
  echo getfilelist $unixtime $trigger failed
  exit 1
fi

$SRT_PRIVATE_CONTEXT/ligo/redoloop_ligo.sh \
strait-ligo-coincidence-artdaq-$unixtime-$trigger
) 2> /dev/stdout | tee $month-$day-$year-$unixtime-$trigger.log
