#!/bin/bash

def="$1"

# totally gratuitous.  Discards stdout of command piped in for first
# 3-4 seconds.  Motivation: "samweb prestage-dataset" says nothing useful
# in that time.
outputafterfirstfewseconds()
{
  t=$(date +%s)
  while read line; do
    tnow=$(date +%s)
    if [ $((tnow - t)) -gt 3 ]; then
      echo "$line"
    fi
  done
}

# Don't know why using $PATH doesn't work...
mycachestate=$SRT_PRIVATE_CONTEXT/NovaGridUtils/bin/cache_state.py
if [ -e $mycachestate ]; then
  CACHESTATE=$mycachestate
else
  CACHESTATE=cache_state.py
fi

# If there's only one file in the set, it says CACHED or NOT
# CACHED.  Otherwise it says "Cached:" and gives a percent.
cachedpercent=$($CACHESTATE -d $def | tee /dev/stderr | \
  awk '/^CACHED$/  {print 100;}\
       /NOT CACHED/{print 0;}\
       /Cached:/   {split($3, n, "("); print n[2]*1;}')


doit()
{
  # Not clear on whether this is necessary, but probably can't hurt
  # Also not clear on, if it is necessary, what the limit should be...
  while true; do
    n=$(ps f -u mstrait | grep 'python2.*prestage-dataset' | grep -v grep  | wc -l)
    if [ $n -gt 3 ]; then
      echo Waiting for $n other prestages to finish
      sleep 1m
    else
      echo Not too many other stages going on.  Going ahead.
      break
    fi
  done

  timeout 10h samweb prestage-dataset --socket-timeout=1800 --defname=$def \
    --parallel 4 2> /dev/stdout | outputafterfirstfewseconds
}

if ! [ $cachedpercent ]; then
  echo Could not see how many files were cached, trying to cache them
  doit
elif [ "$cachedpercent" -lt 100 ]; then
  echo Not all files are cached.  Caching...
  doit
else
  echo All files cached
fi
