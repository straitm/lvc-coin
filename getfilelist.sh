#!/bin/bash -x

echo Using $GWNAME

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

if ! [ $1 ]; then
  echo Specify a trigger stream for second argument
  exit 1
fi

fracsec=$(cut -d. -f 2 -s <<< $gwunixtime)
intsec=$(cut -d. -f 1 <<< $gwunixtime)
rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$intsec).${fracsec}Z
trigger=$1

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
elif [ $trigger == fardet-ddenergy ]; then
  filepattern='fardet.*_ddenergy_.*data.artdaq'
elif [ $trigger == neardet-ddactivity1 ]; then
  # Changed from "DDActivity1" to "ddactivity1" sometime in late 2015, it seems
  filepattern='neardet.*_[Dd][Dd][Aa]ctivity1_.*data.artdaq'
else
  echo unknown trigger \"$trigger\"
  exit 1
fi

defbase=strait-ligo-coincidence-artdaq-$gwunixtime
recodefbase=strait-ligo-coincidence-reco-$gwunixtime
rawdef=$defbase-$trigger
recodef=$recodefbase-$trigger-$gwbase

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
  if samweb describe-definition $def &> /dev/null &&
     ! [ "$(samweb list-files defname: $def)" ]; then
    echo Deleting empty definition
    samweb delete-definition $def
  fi

  if samweb describe-definition $def &> /dev/null; then
    echo SAM definition $def exists for $gwunixtime

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
  # These are the only trigger streams for which we might want to read
  # reco files
  if [ $trigger != fardet-t02 ] &&
     [ $trigger != fardet-ddsnews ] &&
     [ $trigger != neardet-ddactivity1 ]; then
    return 1
  fi

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
    f=$(dirname $outhistdir)/*/*-$trigger/*${base}_*.reco.root
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

  # Some subruns have a start time of zero in the metadata.  In this case,
  # Look at the run start time.  If the run start time isn't there, don't
  # select.  I don't know if this happens. Might there be other problems?
  blocksam
  samweb list-files \
         'Online.SubRunEndTime   > '$((intgwunixtime-550))\
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

echo Running setup_fnal_security.  Will hang if it needs a password.
setup_fnal_security &> /dev/null
echo Ok, it did not hang

# Sleep a little so we can launch a bunch of processes at once without having
# *too* much racing.
if ! [ $REDOFAST ]; then sleep $((RANDOM%16 + 1)); fi

# These get messed up too easily.  For now, regenerate each time.
samweb delete-definition $recodef

if havedef $rawdef; then
  if true; then # XXX
    echo Use artdaq regardless because broken
    def=$rawdef
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

if ! samweb describe-definition $def &> /dev/null; then
  echo Failed to get or make the definition
  exit 1
fi

rm -f $metadir/allfiles.$gwunixtime.$trigger $metadir/selectedfiles.$gwunixtime.$trigger
