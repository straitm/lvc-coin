#!/bin/bash

# Cache the files needed for the 8:30 DDsnews trigger of the given month and day
# in 2019, then run jobs to process them all, running redo loops that look for
# failed jobs and resubmit them.  Try to do the redo loops gently enough that
# we don't fill up /tmp with overlapping submissions too often or cause other
# problems.

md="$1"

unixtime=$(date +%s -d"$md 8:29:01 2019")

if ! $SRT_PRIVATE_CONTEXT/ligo/getfilelist.sh $unixtime; then
  echo getfilelist failed
  exit 1
fi

redoloop_ligo strait-ligo-coincidence-artdaq-$unixtime-neardet-ddsnews &
sleep 2m

redoloop_ligo strait-ligo-coincidence-artdaq-$unixtime-neardet-ddactivity1 &
sleep 2m

redoloop_ligo strait-ligo-coincidence-artdaq-$unixtime-fardet-ddsnews &
sleep 2m

redoloop_ligo strait-ligo-coincidence-artdaq-$unixtime-fardet-ddenergy &
sleep 2m

redoloop_ligo strait-ligo-coincidence-artdaq-$unixtime-fardet-t02 &
sleep 2m
