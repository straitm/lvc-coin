#!/bin/bash

# This gets the time of the GW event, not the time the trigger was issued.

if ! [ $1 ]; then
  echo Give an event name
  exit 1
fi

e=$1

set -o pipefail

url=https://gracedb.ligo.org/apiweb/superevents/$e/files/${e}-1-Preliminary.xml,0

if ! wget --no-check-certificate -q -O - $url | grep ISOTime | cut -d\> -f 2 | cut -d\< -f 1; then
  wget --no-check-certificate $url -O /dev/null
fi
