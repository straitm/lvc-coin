#!/bin/bash

# Unix time stamp of the event we want to search around.
# (How does Online.SubRun*Time handle leap seconds!?)
if [ $1 ]; then
  t=$1
else
  t=1490550400
fi

if ! [ -e allfiles.$t ]; then
  # Fantastically slow
  samweb list-files 'Online.SubRunEndTime > '$((t-500))' and Online.SubRunStartTime < '$((t+500)) > allfiles.$t
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

if [ $2 ] && [ $2 == create-definition ]; then
  def=strait-ligo-coincidence-artdaq-$t
  samweb create-definition $def "$(for f in $(cat selectedfiles.$t | grep -v ^$); do printf "%s %s or " file_name $(basename $f); done | sed 's/or $//')"
  cache_state.py -d $def
  echo consider running: samweb prestage-dataset --defname=$def --parallel 4
fi
