#!/bin/bash

printf '\t\tnd-act\tnd-ligo\tfd-eng\tfd-t02\tfd-ligo\n'

tmpjob=/tmp/joblist.$$

jobsub_q --user mstrait > $tmpjob

noligo=false

no()
{
  printf -- '-\t'
}

for e in \
  S190426c \
  S190513bm \
  S190521g \
  S190602aq \
  S190630ag \
  S190701ah \
  S190706ai \
  \
  GW150914 \
  GW151012 \
  GW151226 \
  GW170104 \
  GW170608 \
  GW170729 \
  GW170809 \
  GW170814 \
  GW170817 \
  GW170818 \
  GW170823 \
  S190412m \
  S190421ar \
  S190425z \
  S190503bf \
  S190510g \
  S190512at \
  S190517h \
  S190519bj \
  S190521r \
  S190707q \
  ; do 
  cd /nova/ana/users/mstrait/ligoresults-$e
  combineall.sh &> /dev/null

  printf $e'\t'
  for trig in neardet-ddactivity1 neardet-ligo fardet-ddenergy fardet-t02 fardet-ligo; do
    if [ -e 2*-${trig}.hadded.root ]; then
      printf 'Y\t'
    elif $noligo && ( [ $trig == fardet-ligo ] || [ $trig == neardet-ligo ] ); then no
    elif [ $e == S190426c  ] && [ $trig == fardet-ligo ]; then no
    elif [ $e == S190513bm ] && [ $trig == fardet-ligo ]; then no
    elif [ $e == S190707q  ] && [ $trig == neardet-ligo ]; then no
    elif [ $e == GW150914   ] && ( [ $trig == fardet-ddenergy ] || [ $trig == fardet-t02 ] ); then no
    elif [ $e == GW151012   ] && ( [ $trig == fardet-ddenergy ] || [ $trig == fardet-t02 ] ); then no
    elif grep -q ${e}'.*\..*'$trig $tmpjob; then
      printf 'Run\t'
    else
      printf 'N\t'
    fi
  done
  printf '\n'

  # Separator between events with long readout and others
  if [ $e == S190706ai ]; then
    printf '\n'
    noligo=true
  fi
done

rm -f $tmpjob
