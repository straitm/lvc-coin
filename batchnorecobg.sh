#!/bin/bash

gweventname=$1

nin=$(find /pnfs/nova/persistent/users/mstrait/ligobg-S190910d/ -name *.reco.root | wc -l)

echo $nin files wanted

find /pnfs/nova/persistent/users/mstrait/ligobg-S190910d/ -name *.reco.root | while read f; do
  rfctime=$(basename $(dirname $f) | cut -d- -f 1-5 | sed -e s/-/:/g -e s/:/-/ -e s/:/-/)
  out=$rfctime-$(basename $f .reco.root).hists.root
  if ! [ -e $out ]; then
    echo Processing $f to $out
    if ! norecobg.sh $gweventname $f; then
      echo Failed to process $f
      rm -fv $out
    fi
  fi
done

nout=$(ls *.hists.root | wc -l)

if [ $nout -eq $nin ]; then
  echo Looks like all files were processed
  hadd -f 109x-fardet-t02.root *t02*hists.root
  hadd -f 229x-neardet-ddactivity1.root *ddactivity1*hists.root
else
  echo Wanted $nin output files, but only got $nout
fi
