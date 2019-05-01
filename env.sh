# Will want a separate directory for each GW's background sample, but for now
# just testing for a single skymap with one.

if ! [ $GWNAME ]; then
  echo export GWNAME=something before calling this
  exit
fi

outhistdir=/pnfs/nova/scratch/users/mstrait/ligobg-$GWNAME/
outhadddir=/nova/ana/users/mstrait/ligobgresults-$GWNAME/

skymap=/pnfs/nova/users/mstrait/ligo/LALInference_skymap-$GWNAME.fits
if [ $GWNAME == S190426c ]; then
  realgweventtime="2019-04-26T15:21:55.3365Z"
else
  echo I do not know when $GWNAME was.  Edit env.sh
  exit 1
fi

if ! echo $PATH | grep -qE "$SRT_PRIVATE_CONTEXT/ligo([:/]|$)"; then
  PATH+=:$SRT_PRIVATE_CONTEXT/ligo
fi

