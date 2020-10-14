#!/bin/bash

if ! [ $1 ]; then
  echo Give an event name
  exit 1
fi

TZ=UTC date +%s -d "$(lvc_event_time.sh $1 | sed 's/T/ /')"
