# Will want a separate directory for each GW's background sample, but for now
# just testing for a single skymap with one.

if ! [ $GWNAME ]; then
  echo export GWNAME=something before calling this
  exit
fi

outhistdir=/pnfs/nova/persistent/users/mstrait/ligobg-$GWNAME/
outhadddir=/nova/ana/users/mstrait/ligobgresults-$GWNAME/

# My naming convention for files I get from GraceDB that don't have the event
# name in them.
skymap=/pnfs/nova/users/mstrait/ligo/LALInference_skymap-$GWNAME.fits

# My naming convention for file I get from GraceDB that are done with bayestar,
# whatever that is.
if ! [ -e $skymap ]; then
  skymap=/pnfs/nova/users/mstrait/ligo/bayestar_skymap-${GWNAME}.fits
fi

# Naming convention for O1 and O2 event catalog
if ! [ -e $skymap ]; then
  skymap=/pnfs/nova/users/mstrait/ligo/${GWNAME}_skymap.fits
fi

# O1 and O2
if   [ $GWNAME == GW150914 ]; then realgweventtime="2015-09-14T09:50:45.4Z"
elif [ $GWNAME == GW151012 ]; then realgweventtime="2015-10-12T09:54:43.4Z"
elif [ $GWNAME == GW151226 ]; then realgweventtime="2015-12-26T03:38:53.65Z"
elif [ $GWNAME == GW170104 ]; then realgweventtime="2017-01-04T10:11:58.6Z"
elif [ $GWNAME == GW170608 ]; then realgweventtime="2017-06-08T02:01:16.49Z"
elif [ $GWNAME == GW170729 ]; then realgweventtime="2017-07-29T18:56:29.3Z"
elif [ $GWNAME == GW170809 ]; then realgweventtime="2017-08-09T08:28:21.8Z"
elif [ $GWNAME == GW170814 ]; then realgweventtime="2017-08-14T10:30:43.53Z"
elif [ $GWNAME == GW170817 ]; then realgweventtime="2017-08-17T12:41:04.4Z"
elif [ $GWNAME == GW170818 ]; then realgweventtime="2017-08-18T02:25:09.1Z"
elif [ $GWNAME == GW170823 ]; then realgweventtime="2017-08-23T13:13:58.5Z"

#03
elif [ $GWNAME == S190412m  ]; then realgweventtime="2019-04-12T05:30:44.1656Z"
elif [ $GWNAME == S190421ar ]; then realgweventtime="2019-04-21T21:38:56.25Z"
elif [ $GWNAME == S190425z  ]; then realgweventtime="2019-04-25T08:08:05.02Z"
elif [ $GWNAME == S190426c  ]; then realgweventtime="2019-04-26T15:21:55.3365Z"
elif [ $GWNAME == S190503bf ]; then realgweventtime="2019-05-03T18:54:04.4126Z"
else
  echo I do not know when $GWNAME was.  Edit env.sh
  if ! [ $PS1 ]; then
    exit 1
  fi
fi

mkdir -p $outhistdir $outhaddir

if ! echo $PATH | grep -qE "$SRT_PRIVATE_CONTEXT/ligo([:/]|$)"; then
  PATH+=:$SRT_PRIVATE_CONTEXT/ligo
fi
