#!/bin/bash

t1=$(date +%s -d $(lvc_event_time.sh $1))
t2=$(date +%s -d $(lvc_trigger_time.sh $1))

echo $t1, $(( (t2 - t1)/60 )).$(( ( (t2 - t1)%60)/6 ))
