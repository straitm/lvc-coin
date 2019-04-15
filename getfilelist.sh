#!/bin/bash

# Unix time stamp of the event we want to search around.
# (How does Online.SubRun*Time handle leap seconds!?)
if ! [ $1 ]; then
  echo Specify a Unix timestamp for first argument
  exit 1
fi

if ! [ $2 ]; then
  echo Specify a trigger stream for second argument
  exit 1
fi

t=$1
trigger=$2

if [ $trigger == neardet-ddsnews ]; then
  filepattern='neardet.*_ddsnews_.*data.artdaq'
elif [ $trigger == fardet-ddsnews ]; then
  filepattern='fardet.*_DDsnews.raw'
elif [ $trigger == neardet-ligo ]; then
  filepattern='neardet.*_ligo_.*data.artdaq'
elif [ $trigger == fardet-ligo ]; then
  filepattern='fardet.*_ligo_.*data.artdaq' # or raw?
elif [ $trigger == fardet-t02 ]; then
  filepattern='fardet.*_t02_.*data.artdaq'
elif [ $trigger == fardet-ddenergy ]; then
  filepattern='fardet.*_ddenergy_.*data.artdaq'
elif [ $trigger == neardet-ddactivity1 ]; then
  filepattern='neardet.*_ddactivity1_.*data.artdaq'
fi

defbase=strait-ligo-coincidence-artdaq-$t

havealldefs()
{
  def=$defbase-$trigger
  if samweb list-definitions | grep -qE "^$def$"; then
    echo SAM definition $def already exists for $t
    return 0
  else
    echo No SAM definition yet for $def
    return 1
  fi
}

nsamlistsrunning()
{
  ps f -u mstrait | grep -v grep | grep samweb.py\ list | wc -l
}

# SAM won't let a user do more than 5 queries at a time.  Self-limit
# to three so that I can also do interactive queries while my scripts run.
blocksam()
{
  while true; do
    n=$(nsamlistsrunning)
    if [ $n -le 2 ]; then break; fi
    echo Waiting for $n sam list processes to finish
    sleep 1m
  done
}

setup_fnal_security &> /dev/null

if havealldefs; then
  echo Have all SAM definitions already.  Doing no queries.
else
  if ! [ -e allfiles.$t.$trigger ]; then
    # I want a 1000 second window, but also add a 50 second buffer on each end.
    # 
    # Also add more time at the end because the Online.SubRun*Time means the
    # times the triggers were issued, not the time of the data.
    # Experimentally, using the 8:30 SNEWS trigger on 2 Feb 2019 as an example,
    # the last SubRunStartTime was 717 seconds after the trigger time, where
    # there is about a 1 minute offset between the trigger time and the data
    # time, so if we could get the trigger as much as 15 minutes after the time
    # it wants, we need 717 - 60 + 15*60 seconds = 1557s, at least. Round up.
    #
    # This also catches cases like how DDEnergy events can be up to 200 seconds
    # late (worse as the events get bigger).
    #
    # Fantastically slow
    echo Asking SAM for a list of files. This typically takes a
    echo substantial fraction of the age of the universe to complete.

    # Some subruns have a start time of zero in the metadata.  In this case,
    # Look at the run start time.  If the run start time isn't there, don't
    # select.  I don't know if this happens. Might there be other problems?
    blocksam
    samweb list-files \
           'Online.SubRunEndTime   > '$((t-550))\
      ' and ( '\
      ' ( Online.SubRunStartTime > 0 and Online.SubRunStartTime < '$((t+2000))\
      ' ) or '\
      ' ( Online.SubRunStartTime = 0 and Online.RunStartTime > 0 and '\
      '                                  Online.RunStartTime < '$((t+2000))\
      ' ) )' \
       > allfiles.$t.$trigger
  fi

  cat allfiles.$t.$trigger | grep -E "$filepattern" | tee selectedfiles.$t.$trigger

  def=$defbase-trigger
  if cat selectedfiles.$t.$trigger | grep -q $filepattern; then
    # Even during the SNEWS trigger, we only get about 40 subruns at the
    # FD in a half hour.  Finding more than that in 1100 seconds (0.3h)
    # means that something is broken.
    if [ $(cat selectedfiles.$t.$trigger|grep $filepattern|wc -l) -gt 99 ]; then
      echo Unreasonable number of files for $trigger.  Skipping.
      continue
    fi

    if samweb list-definitions | grep -qE "^$def$"; then
      continue
    fi

    samweb create-definition $def \
      "$(for f in $(cat selectedfiles.$t.$trigger | grep $filepattern); do
           printf "%s %s or " file_name $(basename $f);
         done | sed 's/ or $//')"

    if ! samweb list-definitions | grep -qE "^$def$"; then
      echo Failed to make a SAM definition
      exit 1
    fi
  fi
fi

def=$defbase-$trigger
if ! samweb list-definitions | grep -qE "^$def$"; then
  continue
fi

$SRT_PRIVATE_CONTEXT/ligo/stage.sh $def

rm -f allfiles.$t.$trigger selectedfiles.$t.$trigger
