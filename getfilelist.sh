#!/bin/bash

nsamlistsrunning()
{
  ps f -u mstrait | grep -v grep | grep samweb.py\ list | wc -l
}

# SAM won't let a user do more than 5 queries at a time.  Self-limit
# to three so that I can also do interactive queries while my scripts run.
blocksam()
{
  if [ $REDOFAST ]; then return; fi
  local try=0
  while true; do
    n=$(nsamlistsrunning)
    if [ $n -le 2 ]; then break; fi
    echo Waiting for $n sam list processes to finish
    sleep $((try + 60 + RANDOM%60))
    let try++
  done
  if [ $try -gt 0 ]; then
    echo Ok, going ahead with SAM query
  fi
}

havedef()
{
  def=$1
  if samweb describe-definition $def &> /dev/null &&
     ! [ "$(samweb list-files defname: $def)" ]; then
    echo Deleting empty definition
    samweb delete-definition $def
  fi

  if samweb describe-definition $def &> /dev/null; then
    echo SAM definition $def exists for $gwunixtime
    return 0
  else
    echo No SAM definition yet for $def
    return 1
  fi
}

makerawdef()
{
  def=$rawdef
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
  echo Asking SAM for a list of files. This typically takes forever.

  # Some subruns have a start time of zero in the metadata.  In this case,
  # Look at the run start time.  If the run start time isn't there, don't
  # select.  I don't know if this happens.  This seems to cover all of the
  # strange cases that happen.
  #
  # Require more than 3 events because that's typically how many SNEWS or
  # GCN (LVC) slow beats there are in files that have no triggered data.
  blocksam
  samweb list-files \
         'Online.TotalEvents > 3 and Online.SubRunEndTime   > '$((intgwunixtime-550))\
    ' and ( '\
    ' ( Online.SubRunStartTime > 0 and Online.SubRunStartTime < '$((intgwunixtime+2000))\
    ' ) or '\
    ' ( Online.SubRunStartTime = 0 and Online.RunStartTime > 0 and '\
    '                                  Online.RunStartTime < '$((intgwunixtime+2000))\
    ' ) )' \
     > $metadir/allfiles.$gwunixtime.$trigger

  echo SAM selected these files:

  cat $metadir/allfiles.$gwunixtime.$trigger | grep -Ei "$filepattern" | sort | \
    tee $metadir/selectedfiles.$gwunixtime.$trigger

  if ! [ "$(cat $metadir/selectedfiles.$gwunixtime.$trigger)" ]; then
    echo No files selected, nothing to do
    if [ $trigger == neardet-t00 ]; then
      echo Making empty spill list file since there are no spills.
      touch $spilldir/spills-$gwunixtime-$rfctime.txt
    fi
    exit 2
  fi

  if cat $metadir/selectedfiles.$gwunixtime.$trigger | grep -qi "$filepattern"; then
    # Even during the SNEWS trigger, we only get about 40 subruns at
    # the FD in a half hour. Finding much more than about three times
    # that number (raw + artdaq + pid) in 1100 seconds (0.3h) means that
    # something is broken.
    if [ $(cat $metadir/selectedfiles.$gwunixtime.$trigger|grep -i "$filepattern"|wc -l) -gt 200 ]; then
      echo Unreasonable number of files for $trigger.
      exit 1
    fi

    # Just in case another script made the definition in the meanwhile?!
    if ! samweb describe-definition $def &> /dev/null; then
      timeout 5m samweb create-definition $def \
        "$(for f in $(grep -i "$filepattern" $metadir/selectedfiles.$gwunixtime.$trigger); do
             printf "file_name %s or " $(basename $f);
           done | sed 's/ or $//')"

      if ! samweb describe-definition $def &> /dev/null; then
        echo Failed to make a SAM definition
        exit 1
      fi
    fi
  fi
}

########################################################################

if ! [ $2 ]; then
  echo Give event name and trigger stream
  exit 1
fi

export GWNAME=$1

echo Using $GWNAME

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

fracsec=$(cut -d. -f 2 -s <<< $gwunixtime)
intsec=$(cut -d. -f 1 <<< $gwunixtime)
rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$intsec).${fracsec}Z
trigger=$2

if [ $trigger == neardet-ddsnews ]; then
  filepattern='neardet.*_DDsnews.raw'
elif [ $trigger == fardet-ddsnews ]; then
  filepattern='fardet.*_DDsnews.raw'
elif [ $trigger == neardet-ligo ]; then
  filepattern='neardet.*_ligo.raw'
elif [ $trigger == fardet-ligo ]; then
  filepattern='fardet.*_ligo.raw'
elif [ $trigger == fardet-ddsn ]; then
  filepattern='fardet.*ddsn'
elif [ $trigger == neardet-ddsn ]; then
  filepattern='neardet.*ddsn'
elif [ $trigger == fardet-t02 ]; then
  filepattern='fardet.*_t02_.*data.artdaq'
elif [ $trigger == neardet-t00 ]; then
  filepattern='neardet.*_t00_.*data.artdaq'
else
  echo unknown trigger \"$trigger\"
  exit 1
fi

defbase=strait-ligo-coincidence-artdaq-$gwunixtime
rawdef=$defbase-$trigger

metadir=/nova/ana/users/mstrait/ligometalog
mkdir -p $metadir

echo Running setup_fnal_security.  Will hang if it needs a password.
setup_fnal_security &> /dev/null
echo Ok, it did not hang

# Sleep a little so we can launch a bunch of processes at once without having
# *too* much racing.
if ! [ $REDOFAST ]; then sleep $((RANDOM%16 + 1)); fi

if ! havedef $rawdef; then
  echo Making raw SAM definition
  makerawdef
  def=$rawdef
fi

if ! samweb describe-definition $def &> /dev/null; then
  echo Failed to get or make the definition
  exit 1
fi

rm -f $metadir/allfiles.$gwunixtime.$trigger $metadir/selectedfiles.$gwunixtime.$trigger
