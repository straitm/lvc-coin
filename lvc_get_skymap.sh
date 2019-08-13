#!/bin/bash

e=$1

if wget -O - https://gracedb.ligo.org/apiweb/superevents/$e/files/LALInference.offline.fits.gz | gunzip > LALInference_skymap-$e.fits; then
  exit 0
elif wget -O - https://gracedb.ligo.org/apiweb/superevents/$e/files/LALInference.fits.gz | gunzip > LALInference_skymap-$e.fits; then
  exit 0
elif wget -O - https://gracedb.ligo.org/apiweb/superevents/$e/files/LALInference1.fits.gz | gunzip > LALInference_skymap-$e.fits; then
  exit 0
else
  rm -f LALInference_skymap-$e.fits
  if wget -O - https://gracedb.ligo.org/apiweb/superevents/$e/files/bayestar.fits.gz | gunzip > bayestar_skymap-$e.fits; then
    exit 0
  else
    rm -f bayestar_skymap-$e.fits
  fi
fi
