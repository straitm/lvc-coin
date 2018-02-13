#!/bin/bash

# Unix time stamp of the event we want to search around.
# (How does Online.SubRun*Time handle leap seconds!?)
if [ $1 ]; then
  t=$1
else
  echo Specify a Unix timestamp
  exit 1
fi

def=strait-ligo-coincidence-artdaq-$t
if samweb list-definitions | grep -qE "^$def$"; then
  echo SAM definition $def already exists for $t
else
  if ! [ -e allfiles.$t ]; then
    # I want a 1000 second window, but also add a 50 second buffer on each
    # end because I have observed, at least with DDEnergy files, that the
    # events are out of order on the scale of 10 seconds and the metadata
    # variables don't seem to know that.
    #
    # Fantastically slow
    echo Asking SAM for a list of files. This typically takes a
    echo stubstantial fraction of the age of the universe to complete.
    samweb list-files \
           'Online.SubRunEndTime   > '$((t-550))\
      ' and Online.SubRunStartTime < '$((t+550)) > allfiles.$t
  fi

  for trigger in \
      'fardet.*_t02_.*data.artdaq' \
      'neardet.*_bnb_.*data.artdaq' \
      'fardet.*_ddenergy_.*data.artdaq' \
      'fardet.*_ddupmu_.*data.artdaq' \
      'fardet.*_ddnumu_.*data.artdaq' \
      'fardet.*_ddfastmono_.*data.artdaq' \
      'fardet.*_ddslowmono_.*data.artdaq' \
      'fardet.*_ddcontained_.*data.artdaq' \
      'neardet.*_ddactivity1_.*data.artdaq' \
      ; do
    cat allfiles.$t | grep -E "$trigger"
    echo
  done | tee selectedfiles.$t

  samweb create-definition $def \
    "$(for f in $(cat selectedfiles.$t | grep -v ^$); do
         printf "%s %s or " file_name $(basename $f);
       done | sed 's/ or $//')"

  if ! samweb list-definitions | grep -qE "^$def$"; then
    echo Failed to make a SAM definition
    exit 1
  fi
fi

cachedpercent=$(cache_state.py -d $def | tee /dev/stderr | \
  awk '/Cached:/{split($1, n, "("); print n[2]*1;}')

if [ "$cachedpercent" -lt 100 ]; then
  echo consider running:
  echo samweb prestage-dataset --defname=$def --parallel 4
fi
