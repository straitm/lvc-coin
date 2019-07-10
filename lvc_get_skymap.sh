#!/bin/bash

e=$1

wget -O - https://gracedb.ligo.org/apiweb/superevents/$e/files/bayestar.fits.gz,0 | gunzip > bayestar_skymap-$e.fits
