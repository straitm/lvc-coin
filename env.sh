if ! [ $GWNAME ]; then
  echo export GWNAME=something before calling this
  exit
fi

gwbase=${GWNAME%%compat}

tag=bg
if [ $REALGWEVENT ]; then
  tag=
fi

outhistdir=/pnfs/nova/persistent/users/mstrait/ligo${tag}-$GWNAME/
outhadddir=/nova/ana/users/mstrait/ligo${tag}results-$GWNAME/
spilldir=/pnfs/nova/users/mstrait/spills

# My naming convention for files I get from GraceDB that don't have the event
# name in them.
skymap=/pnfs/nova/users/mstrait/ligo/LALInference_skymap-$gwbase.fits

# My naming convention for file I get from GraceDB that are done with bayestar,
# whatever that is.
if ! [ -e $skymap ]; then
  skymap=/pnfs/nova/users/mstrait/ligo/bayestar_skymap-${gwbase}.fits
fi

# Naming convention for O1 and O2 event catalog
if ! [ -e $skymap ]; then
  skymap=/pnfs/nova/users/mstrait/ligo/${gwbase}_skymap.fits
fi

# O1 and O2
if   [ $gwbase == GW150914 ]; then realgweventtime="2015-09-14T09:50:45.4Z"
elif [ $gwbase == GW151012 ]; then realgweventtime="2015-10-12T09:54:43.4Z"
elif [ $gwbase == GW151226 ]; then realgweventtime="2015-12-26T03:38:53.65Z"
elif [ $gwbase == GW170104 ]; then realgweventtime="2017-01-04T10:11:58.6Z"
elif [ $gwbase == GW170608 ]; then realgweventtime="2017-06-08T02:01:16.49Z"
elif [ $gwbase == GW170729 ]; then realgweventtime="2017-07-29T18:56:29.3Z"
elif [ $gwbase == GW170809 ]; then realgweventtime="2017-08-09T08:28:21.8Z"
elif [ $gwbase == GW170814 ]; then realgweventtime="2017-08-14T10:30:43.53Z"
elif [ $gwbase == GW170817 ]; then realgweventtime="2017-08-17T12:41:04.4Z"
elif [ $gwbase == GW170818 ]; then realgweventtime="2017-08-18T02:25:09.1Z"
elif [ $gwbase == GW170823 ]; then realgweventtime="2017-08-23T13:13:58.5Z"

#03
elif [ $gwbase == S190412m  ]; then realgweventtime="2019-04-12T05:30:44.165622Z"
elif [ $gwbase == S190421ar ]; then realgweventtime="2019-04-21T21:38:56.250977Z"
elif [ $gwbase == S190425z  ]; then realgweventtime="2019-04-25T08:18:05.017147Z"
elif [ $gwbase == S190426c  ]; then realgweventtime="2019-04-26T15:21:55.336540Z"
elif [ $gwbase == S190503bf ]; then realgweventtime="2019-05-03T18:54:04.294490Z"
elif [ $gwbase == S190510g  ]; then realgweventtime="2019-05-10T02:59:39.291636Z"
elif [ $gwbase == S190512at ]; then realgweventtime="2019-05-12T18:07:14.422363Z"
elif [ $gwbase == S190513bm ]; then realgweventtime="2019-05-13T20:54:28.747089Z"
elif [ $gwbase == S190517h  ]; then realgweventtime="2019-05-17T05:51:01.830582Z"
elif [ $gwbase == S190519bj ]; then realgweventtime="2019-05-19T15:35:44.397949Z"
elif [ $gwbase == S190521g  ]; then realgweventtime="2019-05-21T03:02:29.447266Z"
elif [ $gwbase == S190521r  ]; then realgweventtime="2019-05-21T07:43:59.463379Z"
else
  echo I do not know when $gwbase was.  Edit env.sh
  if ! [ $PS1 ]; then
    exit 1
  fi
fi

# For comparison with older tests
if [ $GWCOMPAT ]; then
  realgweventtime="2019-03-24T13:29:01Z"
  GWNAME=${gwbase}compat
fi

mkdir -p $outhistdir $outhaddir

if ! echo $PATH | grep -qE "$SRT_PRIVATE_CONTEXT/ligo([:/]|$)"; then
  PATH+=:$SRT_PRIVATE_CONTEXT/ligo
fi
