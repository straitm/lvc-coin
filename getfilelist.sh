#!/bin/bash

# Prints a random 32 bit unsigned number.  Not cryptographically secure.
# Probably not even good enough for serious Monte Carlo, but good enough
# to sample some data.
rand32()
{
  echo $(( $RANDOM * 0x10000 + $RANDOM))
}

# Pick a random Unix time stamp (or, really, whatever number) between $1
# and $2. 32-bit *unsigned* numbers.  Fails after Jan 2106.
randtime()
{
  start=$1
  end=$2

  # Ultra-stupid retry loop
  # Did you know you could write while loops in bash like this?
  while
    r=$(rand32)
    [ $r -lt $start ] || [ $r -gt $end ]
  do
    :
  done
  echo $r
}

# Unix time stamp of the event we want to search around.
# (How does Online.SubRun*Time handle leap seconds!?)
if [ $1 ]; then
  if [ "$1" == random ]; then
    start=1359698400
    end=$(date +%s)
    t=$(randtime $start $end)
  else
    t=$1
  fi
else
  echo Specify a Unix timestamp
  exit 1
fi

# Change "ddsnews" to "ligo" when ligo files are available
filepatterns=(
  'neardet.*_ddsnews_.*data.artdaq'
  'fardet.*_DDsnews.raw'
  'fardet.*_t02_.*data.artdaq'
  'fardet.*_ddenergy_.*data.artdaq'
  'neardet.*_ddactivity1_.*data.artdaq'
)

# Change "ddsnews" to "ligo" when ligo files are available
triggers=(
  'neardet-ddsnews'
  'fardet-ddsnews'
  'fardet-t02'
  'fardet-ddenergy'
  'neardet-ddactivity1'
)

defbase=strait-ligo-coincidence-artdaq-$t

havealldefs()
{
  need=0
  for i in {0..4}; do
    def=$defbase-${triggers[i]}
    if samweb list-definitions | grep -qE "^$def$"; then
      echo SAM definition $def already exists for $t
    else
      need=1
      echo No SAM definition yet for $def
    fi
  done
  return $need
}

setup_fnal_security &> /dev/null

if havealldefs; then
  echo Have all SAM definitions already.  Doing no queries.
else
  if ! [ -e allfiles.$t ]; then
    # I want a 1000 second window, but also add a 50 second buffer on each
    # end because I have observed, at least with DDEnergy files, that the
    # events are out of order on the scale of 10 seconds and the metadata
    # variables don't seem to know that.
    #
    # Also add more time at the end because the Online.SubRun*Time means
    # the times the triggers were issued, not the time of the data.  Experimentally,
    # using the 8:30 SNEWS trigger on 2 Feb 2019 as an example, the last
    # SubRunStartTime was 717 seconds after the trigger time, where there is
    # about a 1 minute offset between the trigger time and the data time, so
    # if we could get the trigger as much as 15 minutes after the time it
    # wants, we need 717 - 60 + 15*60 seconds = 1557s, at least. Round up.
    #
    # Fantastically slow
    echo Asking SAM for a list of files. This typically takes a
    echo substantial fraction of the age of the universe to complete.

    # Some subruns have a start time of zero in the metadata.  In this case,
    # Look at the run start time.  If the run start time isn't there, don't
    # select.  I don't know if this happens. Might there be other problems?
    samweb list-files \
           'Online.SubRunEndTime   > '$((t-550))\
      ' and ( '\
      ' ( Online.SubRunStartTime > 0 and Online.SubRunStartTime < '$((t+2000))\
      ' ) or '\
      ' ( Online.SubRunStartTime = 0 and Online.RunStartTime > 0 and '\
      '                                  Online.RunStartTime < '$((t+2000))\
      ' ) )' \
       > allfiles.$t
  fi

  for i in {0..4}; do
    cat allfiles.$t | grep -E "${filepatterns[i]}"
    echo
  done | tee selectedfiles.$t

  for i in {0..4}; do
    def=$defbase-${triggers[i]}
    if cat selectedfiles.$t | grep -q ${filepatterns[i]}; then
      # Even during the SNEWS trigger, we only get about 40 subruns at the
      # FD in a half hour.  Finding more than that in 1100 seconds (0.3h)
      # means that something is broken.
      if [ $(cat selectedfiles.$t|grep ${filepatterns[i]}|wc -l) -gt 99 ]; then
        echo Unreasonable number of files for ${triggers[i]}.  Skipping.
        continue
      fi

      if samweb list-definitions | grep -qE "^$def$"; then
        continue
      fi

      samweb create-definition $def \
        "$(for f in $(cat selectedfiles.$t | grep ${filepatterns[i]}); do
             printf "%s %s or " file_name $(basename $f);
           done | sed 's/ or $//')"

      if ! samweb list-definitions | grep -qE "^$def$"; then
        echo Failed to make a SAM definition
        exit 1
      fi
    fi
  done
fi

# totally gratuitous.  Discards stdout of command piped in for first
# 1-2 seconds.  Motivation: "samweb prestage-dataset" says nothing useful
# in that time.
outputafterfirstsecond()
{
  t=$(date +%s)
  while read line; do
    tnow=$(date +%s)
    if [ $((tnow - t)) -gt 1 ]; then
      echo "$line"
    fi  
  done
}

for i in {0..4}; do
  def=$defbase-${triggers[i]}
  if ! samweb list-definitions | grep -qE "^$def$"; then
    continue
  fi

  # If there's only one file in the set, it says CACHED or NOT
  # CACHED.  Otherwise it says "Cached:" and gives a percent.
  cachedpercent=$(cache_state.py -d $def | tee /dev/stderr | \
    awk '/^CACHED$/  {print 100;}\
         /NOT CACHED/{print 0;}\
         /Cached:/   {split($3, n, "("); print n[2]*1;}')

  if ! [ $cachedpercent ]; then
    echo Could not see how many files were cached
    continue
  fi

  if [ "$cachedpercent" -lt 100 ]; then
    echo Not all files are cached.  Caching...
    samweb prestage-dataset --defname=$def --parallel 4 | outputafterfirstsecond
  fi
done

rm -f allfiles.$t selectedfiles.$t
