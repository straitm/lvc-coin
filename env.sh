if ! [ $GWNAME ]; then
  echo export GWNAME=something before calling this
  exit
fi

gwbase=${GWNAME%%compat}

tag=bg
if [ $SIDEBAND ]; then
  tag=sideband
elif [ $REALGWEVENT ]; then
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

# From PHYS. REV. D 100, 023007 (2019)
# No skymap available?
elif [ $gwbase == GW151216  ]; then realgweventtime="2015-12-16T09:24:16.165Z"

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
elif [ $gwbase == S190602aq ]; then realgweventtime="2019-06-02T17:59:27.089355Z"
elif [ $gwbase == S190630ag ]; then realgweventtime="2019-06-30T18:52:05.179550Z"
elif [ $gwbase == S190701ah ]; then realgweventtime="2019-07-01T20:33:06.577637Z"
elif [ $gwbase == S190706ai ]; then realgweventtime="2019-07-06T22:26:41.344727Z"
elif [ $gwbase == S190707q  ]; then realgweventtime="2019-07-07T09:33:26.181226Z"

elif [ $gwbase == S190718y  ]; then realgweventtime="2019-07-18T14:35:12.067865Z"
elif [ $gwbase == S190720a  ]; then realgweventtime="2019-07-20T00:08:36.747070Z"
elif [ $gwbase == S190727h  ]; then realgweventtime="2019-07-27T06:03:33.985887Z"
elif [ $gwbase == S190728q  ]; then realgweventtime="2019-07-28T06:45:10.546797Z"
elif [ $gwbase == S190814bv ]; then realgweventtime="2019-08-14T21:10:39.013334Z"
elif [ $gwbase == S190828j  ]; then realgweventtime="2019-08-28T06:34:05.756472Z"
elif [ $gwbase == S190828l  ]; then realgweventtime="2019-08-28T06:55:09.886557Z"
elif [ $gwbase == S190901ap ]; then realgweventtime="2019-09-01T23:31:01.837767Z"
elif [ $gwbase == S190910d  ]; then realgweventtime="2019-09-10T01:26:19.242676Z"
elif [ $gwbase == S190910h  ]; then realgweventtime="2019-09-10T08:29:58.544448Z"
elif [ $gwbase == S190915ak ]; then realgweventtime="2019-09-15T23:57:02.690891Z"
elif [ $gwbase == S190923y  ]; then realgweventtime="2019-09-23T12:55:59.645508Z"
elif [ $gwbase == S190924h  ]; then realgweventtime="2019-09-24T02:18:46.846654Z"
elif [ $gwbase == S190930s  ]; then realgweventtime="2019-09-30T13:35:41.246810Z"
elif [ $gwbase == S190930t  ]; then realgweventtime="2019-09-30T14:34:07.685342Z"

elif [ $gwbase == S200316bj ]; then realgweventtime="2020-03-16T21:57:56.157221Z"

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

mkdir -p $outhistdir $outhadddir

if ! echo $PATH | grep -qE "$SRT_PRIVATE_CONTEXT/ligo([:/]|$)"; then
  PATH+=:$SRT_PRIVATE_CONTEXT/ligo
fi
