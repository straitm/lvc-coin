#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

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

unixtime=$1
fracsec=$(cut -d. -f 2 -s <<< $unixtime)
intsec=$(cut -d. -f 1 <<< $unixtime)
rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$intsec).${fracsec}Z
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
elif [ $trigger == neardet-t00 ]; then
  filepattern='neardet.*_t00_.*data.artdaq'
elif [ $trigger == fardet-ddenergy ]; then
  filepattern='fardet.*_ddenergy_.*data.artdaq'
elif [ $trigger == neardet-ddactivity1 ]; then
  filepattern='neardet.*_ddactivity1_.*data.artdaq'
else
  echo unknown trigger \"$trigger\"
  exit 1
fi

defbase=strait-ligo-coincidence-artdaq-$unixtime
recodefbase=strait-ligo-coincidence-reco-$unixtime
rawdef=$defbase-$trigger
recodef=$recodefbase-$trigger

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
  tmplist=/tmp/samlist.$$

  def=$1
  blocksam
  if samweb list-definitions | grep -qE "^$def$" &&
     ! [ "$(samweb list-files defname: $def)" ]; then
    echo Deleting empty definition
    samweb delete-definition $def
  fi

  blocksam
  if samweb list-definitions | grep -qE "^$def$"; then
    echo SAM definition $def exists for $unixtime

    samweb list-files defname: $def > $tmplist
    for f in $(cat $tmplist); do
      if samweb locate-file $f | grep -q persistent && 
         ! [ -e $(samweb locate-file $f | sed s/dcache://)/$f ]; then
        echo $f in reco def does not exist.  Removing definition.
        samweb delete-definition $recodef
        rm -f $tmplist
        return 1
      fi
    done
    rm -f $tmplist
    return 0
  else
    echo No SAM definition yet for $def
    return 1
  fi
}

makerecodef()
{
  tmplist=/tmp/tmplist.$$
  tmprecolist=/tmp/recolist.$$
  rm -f $tmprecolist

  if ! havedef $rawdef > /dev/null; then
    echo Will make the raw definition before making reco def
    makerawdef
  fi

  blocksam
  samweb list-files defname: $rawdef > $tmplist
  echo Looking for $(cat $tmplist | wc -l) reco file'(s)'...
  for raw in $(cat $tmplist); do
    base=$(printf $raw | cut -d_ -f 1-4 | cut -d. -f 1 | sed s/DDsnews/ddsnews/)
    f=$outhistdir/../*/*/*${base}_*.reco.root 
    echo Looking for "$f"
    if ls $f &> /dev/null; then
      found=0
      for one in $(ls $f); do
        rfctimeonfile=$(basename $(dirname $one) | cut -dZ -f 1)Z
        rfctimeonfile=${rfctimeonfile//-/:} #repairing...
        rfctimeonfile=${rfctimeonfile/:/-}
        rfctimeonfile=${rfctimeonfile/:/-}
        if [ $rfctime == $rfctimeonfile ]; then
          ls $one >> $tmprecolist
          echo Found $one
          found=1
          break
        else
          echo Incompatible timestamp on already-filtered reco file $one
        fi
      done
      if [ $found -ne 1 ]; then
        echo No acceptable reco file found
        rm -f $tmplist $tmprecolist
        return 1
      fi
    else
      echo No reco file found
      rm -f $tmplist $tmprecolist
      return 1
    fi
  done

  if [ -e $tmprecolist ]; then
    for f in $(cat $tmprecolist); do # ???
      samweb retire-file $(basename $f) &> /dev/null
    done
    sam_add_dataset -n $recodef -f $tmprecolist
  else
    return 1
  fi

  rm -f $tmplist $tmprecolist
}

metadir=/nova/ana/users/mstrait/ligometalog
mkdir -p $metadir

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

  # Remove fractional part for integer arithmetic below.  Sub-second
  # precision isn't needed since we leave large buffers on either side.
  intunixtime=${unixtime%.*}

  # Some subruns have a start time of zero in the metadata.  In this case,
  # Look at the run start time.  If the run start time isn't there, don't
  # select.  I don't know if this happens. Might there be other problems?
  blocksam
  samweb list-files \
         'Online.SubRunEndTime   > '$((intunixtime-550))\
    ' and ( '\
    ' ( Online.SubRunStartTime > 0 and Online.SubRunStartTime < '$((intunixtime+2000))\
    ' ) or '\
    ' ( Online.SubRunStartTime = 0 and Online.RunStartTime > 0 and '\
    '                                  Online.RunStartTime < '$((intunixtime+2000))\
    ' ) )' \
     > $metadir/allfiles.$unixtime.$trigger

  echo SAM selected these files:

  cat $metadir/allfiles.$unixtime.$trigger | grep -E "$filepattern" | sort | \
    tee $metadir/selectedfiles.$unixtime.$trigger

  if ! [ "$(cat $metadir/selectedfiles.$unixtime.$trigger)" ]; then
    echo No files selected, nothing to do
    if [ $trigger == neardet-t00 ]; then
      touch $spilldir/spills-$unixtime-$rfctime.txt
    fi
    exit 2
  fi

  if cat $metadir/selectedfiles.$unixtime.$trigger | grep -q "$filepattern"; then
    # Even during the SNEWS trigger, we only get about 40 subruns at the
    # FD in a half hour.  Finding more than that in 1100 seconds (0.3h)
    # means that something is broken.
    if [ $(cat $metadir/selectedfiles.$unixtime.$trigger|grep "$filepattern"|wc -l) -gt 99 ]; then
      echo Unreasonable number of files for $trigger.
      exit 1
    fi

    # Just in case another script made the definition in the meanwhile?!
    blocksam
    if ! samweb list-definitions | grep -qE "^$def$"; then
      samweb create-definition $def \
        "$(for f in $(cat $metadir/selectedfiles.$unixtime.$trigger | grep "$filepattern"); do
             printf "%s %s or " file_name $(basename $f);
           done | sed 's/ or $//')"

      if ! samweb list-definitions | grep -qE "^$def$"; then
        echo Failed to make a SAM definition
        exit 1
      fi
    fi
  fi
}

setup_fnal_security &> /dev/null

# Sleep a little so we can launch a bunch of processes at once without having
# *too* much racing.
if ! [ $REDOFAST ]; then sleep $((RANDOM%16 + 1)); fi

if havedef $rawdef; then
  if havedef $recodef; then
    echo Have raw and reco SAM definitions already.  Doing no queries.
    def=$recodef
  elif makerecodef; then
    echo Made reco SAM definition
    def=$recodef
  else
    echo Reco files not available. Have raw SAM definition. Doing no queries.
    def=$rawdef
  fi
else
  echo Making raw SAM definition
  makerawdef
  def=$rawdef
fi

blocksam
if ! samweb list-definitions | grep -qE "^$def$"; then
  echo Failed to get or make the definition
  exit 1
fi

rm -f $metadir/allfiles.$unixtime.$trigger $metadir/selectedfiles.$unixtime.$trigger
