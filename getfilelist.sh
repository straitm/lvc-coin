#!/bin/bash

# Unix time stamp of the event we want to search around.
# (How does Online.SubRun*Time handle leap seconds!?)
if [ $1 ]; then
  t=$1
else
  echo Specify a Unix timestamp
  exit 1
fi

filepatterns=(
  'fardet.*_t02_.*data.artdaq'
  'neardet.*_bnb_.*data.artdaq'
  'fardet.*_ddenergy_.*data.artdaq'
  'fardet.*_ddupmu_.*data.artdaq'
  'fardet.*_ddnumu_.*data.artdaq'
  'fardet.*_ddfastmono_.*data.artdaq'
  'fardet.*_ddslowmono_.*data.artdaq'
  'fardet.*_ddcontained_.*data.artdaq'
  'neardet.*_ddactivity1_.*data.artdaq'
)

triggers=(
  'fardet-t02'
  'neardet-bnb'
  'fardet-ddenergy'
  'fardet-ddupmu'
  'fardet-ddnumu'
  'fardet-ddfastmono'
  'fardet-ddslowmono'
  'fardet-ddcontained'
  'neardet-ddactivity1'
)

defbase=strait-ligo-coincidence-artdaq-$t

havealldefs()
{
  need=0
  for i in {0..8}; do
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

if havealldefs; then
  echo Have all SAM definitions already
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

    # Some subruns have a start time of zero in the metadata.  In this case,
    # Look at the run start time.  If the run start time isn't there, don't
    # select.  I don't know if this happens. Might there be other problems?
    samweb list-files \
           'Online.SubRunEndTime   > '$((t-550))\
      ' and ( '\
      ' ( Online.SubRunStartTime > 0 and Online.SubRunStartTime < '$((t+550))\
      ' ) or '\
      ' ( Online.SubRunStartTime = 0 and Online.RunStartTime > 0 and '\
      '                                  Online.RunStartTime < '$((t+550))\
      ' ) )' \
       > allfiles.$t
  fi

  for i in {0..8}; do
    cat allfiles.$t | grep -E "${filepatterns[i]}"
    echo
  done | tee selectedfiles.$t

  for i in {0..8}; do
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

for i in {0..8}; do
  def=$defbase-${triggers[i]}
  if ! samweb list-definitions | grep -qE "^$def$"; then
    continue
  fi
  cachedpercent=$(cache_state.py -d $def | tee /dev/stderr | \
    awk '/Cached:/{split($3, n, "("); print n[2]*1;}')

  if [ "$cachedpercent" -lt 100 ]; then
    echo Not all files are cached.  Caching...
    samweb prestage-dataset --defname=$def --parallel 4
  fi
done
