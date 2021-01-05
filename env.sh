if ! [ $GWNAME ]; then
  echo export GWNAME=something before calling this
else
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
  elif [ $gwbase == S190412m  ] || [ $gwbase == GW190412 ]; then realgweventtime="2019-04-12T05:30:44.165622Z"
  elif [ $gwbase == S190421ar ] || [ $gwbase == GW190421_213856 ]; then realgweventtime="2019-04-21T21:38:56.250977Z"
  elif [ $gwbase == S190425z  ] || [ $gwbase == GW190425 ]; then realgweventtime="2019-04-25T08:18:05.017147Z"
  elif [ $gwbase == S190426c  ] || [ $gwbase == GW190426_152155 ]; then realgweventtime="2019-04-26T15:21:55.336540Z"
  elif [ $gwbase == S190503bf ] || [ $gwbase == GW190503_185404 ]; then realgweventtime="2019-05-03T18:54:04.294490Z"
  elif [ $gwbase == S190510g  ]; then realgweventtime="2019-05-10T02:59:39.291636Z"
  elif [ $gwbase == S190512at ] || [ $gwbase == GW190512_180714 ]; then realgweventtime="2019-05-12T18:07:14.422363Z"
  elif [ $gwbase == S190513bm ] || [ $gwbase == GW190513_205428 ]; then realgweventtime="2019-05-13T20:54:28.747089Z"
  elif [ $gwbase == S190517h  ] || [ $gwbase == GW190517_055101 ]; then realgweventtime="2019-05-17T05:51:01.830582Z"
  elif [ $gwbase == S190519bj ] || [ $gwbase == GW190519_153544 ]; then realgweventtime="2019-05-19T15:35:44.397949Z"
  elif [ $gwbase == S190521g  ] || [ $gwbase == GW190521 ]; then realgweventtime="2019-05-21T03:02:29.447266Z"
  elif [ $gwbase == S190521r  ] || [ $gwbase == GW190521_074359 ]; then realgweventtime="2019-05-21T07:43:59.463379Z"
  elif [ $gwbase == S190602aq ] || [ $gwbase == GW190602_175927 ]; then realgweventtime="2019-06-02T17:59:27.089355Z"
  elif [ $gwbase == S190630ag ] || [ $gwbase == GW190630_185205 ]; then realgweventtime="2019-06-30T18:52:05.179550Z"
  elif [ $gwbase == S190701ah ] || [ $gwbase == GW190701_203306 ]; then realgweventtime="2019-07-01T20:33:06.577637Z"
  elif [ $gwbase == S190706ai ] || [ $gwbase == GW190706_222641 ]; then realgweventtime="2019-07-06T22:26:41.344727Z"
  elif [ $gwbase == S190707q  ] || [ $gwbase == GW190707_093326 ]; then realgweventtime="2019-07-07T09:33:26.181226Z"

  # Still O3, my new analysis
  elif [ $gwbase == S190718y  ]; then realgweventtime="2019-07-18T14:35:12.067865Z"
  elif [ $gwbase == S190720a  ] || [ $gwbase == GW190720_000836 ]; then realgweventtime="2019-07-20T00:08:36.747070Z"
  elif [ $gwbase == S190727h  ] || [ $gwbase == GW190727_060333 ]; then realgweventtime="2019-07-27T06:03:33.985887Z"
  elif [ $gwbase == S190728q  ] || [ $gwbase == GW190728_064510 ]; then realgweventtime="2019-07-28T06:45:10.546797Z"
  elif [ $gwbase == S190814bv ] || [ $gwbase == GW190814 ]; then realgweventtime="2019-08-14T21:10:39.013334Z"
  elif [ $gwbase == S190828j  ] || [ $gwbase == GW190828_063405 ]; then realgweventtime="2019-08-28T06:34:05.756472Z"
  elif [ $gwbase == S190828l  ] || [ $gwbase == GW190828_065509 ]; then realgweventtime="2019-08-28T06:55:09.886557Z"
  elif [ $gwbase == S190901ap ]; then realgweventtime="2019-09-01T23:31:01.837767Z"
  elif [ $gwbase == S190910d  ]; then realgweventtime="2019-09-10T01:26:19.242676Z"
  elif [ $gwbase == S190910h  ]; then realgweventtime="2019-09-10T08:29:58.544448Z"
  elif [ $gwbase == S190915ak ] || [ $gwbase == GW190915_235702 ]; then realgweventtime="2019-09-15T23:57:02.690891Z"
  elif [ $gwbase == S190923y  ]; then realgweventtime="2019-09-23T12:55:59.645508Z"
  elif [ $gwbase == S190924h  ] || [ $gwbase == GW190924_021846 ]; then realgweventtime="2019-09-24T02:18:46.846654Z"
  elif [ $gwbase == S190930s  ] || [ $gwbase == GW190930_133541 ]; then realgweventtime="2019-09-30T13:35:41.246810Z"
  elif [ $gwbase == S190930t  ]; then realgweventtime="2019-09-30T14:34:07.685342Z"
  elif [ $gwbase == S191105e  ]; then realgweventtime="2019-11-05T14:35:21.933105Z"
  elif [ $gwbase == S191109d  ]; then realgweventtime="2019-11-09T01:07:17.220703Z"
  elif [ $gwbase == S191129u  ]; then realgweventtime="2019-11-29T13:40:29.197372Z"
  elif [ $gwbase == S191204r  ]; then realgweventtime="2019-12-04T17:15:26.091822Z"
  elif [ $gwbase == S191205ah ]; then realgweventtime="2019-12-05T21:52:08.568738Z"
  elif [ $gwbase == S191213g  ]; then realgweventtime="2019-12-13T04:34:08.142224Z"
  elif [ $gwbase == S191215w  ]; then realgweventtime="2019-12-15T22:30:52.333152Z"
  elif [ $gwbase == S191216ap ]; then realgweventtime="2019-12-16T21:33:38.472999Z"
  elif [ $gwbase == S191222n  ]; then realgweventtime="2019-12-22T03:35:37.119478Z"
  elif [ $gwbase == S200105ae ]; then realgweventtime="2020-01-05T16:24:26.057208Z"
  elif [ $gwbase == S200112r  ]; then realgweventtime="2020-01-12T15:58:38.093931Z"
  elif [ $gwbase == S200114f  ]; then realgweventtime="2020-01-14T02:08:18.239300Z"
  elif [ $gwbase == S200115j  ]; then realgweventtime="2020-01-15T04:23:09.752869Z"
  elif [ $gwbase == S200128d  ]; then realgweventtime="2020-01-28T02:20:11.903320Z"
  elif [ $gwbase == S200129m  ]; then realgweventtime="2020-01-29T06:54:58.435104Z"
  elif [ $gwbase == S200208q  ]; then realgweventtime="2020-02-08T13:01:17.991118Z"
  elif [ $gwbase == S200213t  ]; then realgweventtime="2020-02-13T04:10:40.327981Z"
  elif [ $gwbase == S200219ac ]; then realgweventtime="2020-02-19T09:44:15.195312Z"
  elif [ $gwbase == S200224ca ]; then realgweventtime="2020-02-24T22:22:34.378418Z"
  elif [ $gwbase == S200225q  ]; then realgweventtime="2020-02-25T06:04:21.396484Z"
  elif [ $gwbase == S200302c  ]; then realgweventtime="2020-03-02T01:58:11.519119Z"
  elif [ $gwbase == S200311bg ]; then realgweventtime="2020-03-11T11:58:53.397788Z"
  elif [ $gwbase == S200316bj ]; then realgweventtime="2020-03-16T21:57:56.157221Z"

  # #!/bin/bash
  # 
  # in=$1
  # 
  # year=20${in:2:2}
  # month=${in:4:2}
  # day=${in:6:2}
  # hour=${in:9:2}
  # minute=${in:11:2}
  # second=${in:13:2}
  # 
  # echo 'elif [ '\$gwbase' == '$in' ]; then realgweventtime="'$year-$month-${day}T$hour:$minute:${second}.0Z\"

  elif [ $gwbase == GW190413_052954 ]; then realgweventtime="2019-04-13T05:29:54.0Z"
  elif [ $gwbase == GW190413_134308 ]; then realgweventtime="2019-04-13T13:43:08.0Z"
  elif [ $gwbase == GW190424_180648 ]; then realgweventtime="2019-04-24T18:06:48.0Z"
  elif [ $gwbase == GW190514_065416 ]; then realgweventtime="2019-05-14T06:54:16.0Z"
  elif [ $gwbase == GW190527_092055 ]; then realgweventtime="2019-05-27T09:20:55.0Z"
  elif [ $gwbase == GW190620_030421 ]; then realgweventtime="2019-06-20T03:04:21.0Z"
  elif [ $gwbase == GW190708_232457 ]; then realgweventtime="2019-07-08T23:24:57.0Z"
  elif [ $gwbase == GW190719_215514 ]; then realgweventtime="2019-07-19T21:55:14.0Z"
  elif [ $gwbase == GW190731_140936 ]; then realgweventtime="2019-07-31T14:09:36.0Z"
  elif [ $gwbase == GW190803_022701 ]; then realgweventtime="2019-08-03T02:27:01.0Z"
  elif [ $gwbase == GW190910_112807 ]; then realgweventtime="2019-09-10T11:28:07.0Z"
  elif [ $gwbase == GW190929_012149 ]; then realgweventtime="2019-09-29T01:21:49.0Z"


  # SNEWS test trigger 8:30 Fermilab time = 13:29:01.35 during daylight savings time
  #                                         14:29:01.35 during standard time
  # That is, recently, since we changed the offset, whenever that was.

  elif [ $gwbase == snews20191201 ]; then realgweventtime="2019-12-01T14:29:01.350000Z"

  elif [ $gwbase == snews20200101 ]; then realgweventtime="2020-01-01T14:29:01.350000Z"
  elif [ $gwbase == snews20200115 ]; then realgweventtime="2020-01-15T14:29:01.350000Z"
  elif [ $gwbase == snews20200201 ]; then realgweventtime="2020-02-01T14:29:01.350000Z"
  elif [ $gwbase == snews20200215 ]; then realgweventtime="2020-02-15T14:29:01.350000Z"
  elif [ $gwbase == snews20200301 ]; then realgweventtime="2020-03-01T14:29:01.350000Z"

  # DST began on March 8, 2020
  elif [ $gwbase == snews20200315 ]; then realgweventtime="2020-03-15T14:29:01.350000Z"

  elif [ $gwbase == snews20200401 ]; then realgweventtime="2020-04-01T13:29:01.350000Z"
  elif [ $gwbase == snews20200415 ]; then realgweventtime="2020-04-15T13:29:01.350000Z"
  elif [ $gwbase == snews20200501 ]; then realgweventtime="2020-05-01T13:29:01.350000Z"
  elif [ $gwbase == snews20200515 ]; then realgweventtime="2020-05-15T13:29:01.350000Z"
  elif [ $gwbase == snews20200601 ]; then realgweventtime="2020-06-01T13:29:01.350000Z"
  elif [ $gwbase == snews20200615 ]; then realgweventtime="2020-06-15T13:29:01.350000Z"

  # SNEWS test trigger used for my training test
  elif [ $gwbase == snews20200623 ]; then realgweventtime="2020-06-23T13:29:01.350000Z"

  elif [ $gwbase == snews20200701 ]; then realgweventtime="2020-07-01T13:29:01.350000Z"

  # SNEWS test triggers used for my training test
  elif [ $gwbase == snews20200707 ]; then realgweventtime="2020-07-07T13:29:01.350000Z"
  elif [ $gwbase == snews20200709 ]; then realgweventtime="2020-07-09T13:29:01.350000Z"
  elif [ $gwbase == snews20200710 ]; then realgweventtime="2020-07-10T13:29:01.350000Z"

  elif [ $gwbase == snews20200715 ]; then realgweventtime="2020-07-15T13:29:01.350000Z"

  elif [ $gwbase == snews20200801 ]; then realgweventtime="2020-08-01T13:29:01.350000Z"
  elif [ $gwbase == snews20200815 ]; then realgweventtime="2020-08-15T13:29:01.350000Z"
  elif [ $gwbase == snews20200901 ]; then realgweventtime="2020-09-01T13:29:01.350000Z"
  elif [ $gwbase == snews20200915 ]; then realgweventtime="2020-09-15T13:29:01.350000Z"
  elif [ $gwbase == snews20201001 ]; then realgweventtime="2020-10-01T13:29:01.350000Z"

  # A recent run when I took it
  elif [ $gwbase == snews20201012 ]; then realgweventtime="2020-10-12T13:29:01.350000Z"

  elif [ $gwbase == snews20201015 ]; then realgweventtime="2020-10-15T13:29:01.350000Z"

  # DST ended on Nov 1, 2020

  elif [ $gwbase == snews20201101 ]; then realgweventtime="2020-11-01T14:29:01.350000Z"
  elif [ $gwbase == snews20201102 ]; then realgweventtime="2020-11-02T14:29:01.350000Z"
  elif [ $gwbase == snews20201103 ]; then realgweventtime="2020-11-03T14:29:01.350000Z"

  elif [ $gwbase == snews20201104 ]; then realgweventtime="2020-11-04T14:29:01.350000Z"
  elif [ $gwbase == snews20201105 ]; then realgweventtime="2020-11-05T14:29:01.350000Z"
  elif [ $gwbase == snews20201106 ]; then realgweventtime="2020-11-06T14:29:01.350000Z"
  elif [ $gwbase == snews20201107 ]; then realgweventtime="2020-11-07T14:29:01.350000Z"
  elif [ $gwbase == snews20201108 ]; then realgweventtime="2020-11-08T14:29:01.350000Z"
  elif [ $gwbase == snews20201109 ]; then realgweventtime="2020-11-09T14:29:01.350000Z"

  elif [ $gwbase == snews20201110 ]; then realgweventtime="2020-11-10T14:29:01.350000Z"
  elif [ $gwbase == snews20201111 ]; then realgweventtime="2020-11-11T14:29:01.350000Z"
  elif [ $gwbase == snews20201112 ]; then realgweventtime="2020-11-12T14:29:01.350000Z"
  elif [ $gwbase == snews20201113 ]; then realgweventtime="2020-11-13T14:29:01.350000Z"
  elif [ $gwbase == snews20201114 ]; then realgweventtime="2020-11-14T14:29:01.350000Z"

  # Used for overlaying MC
  elif [ $gwbase == snews20201115 ]; then realgweventtime="2020-11-15T14:29:01.350000Z"

  elif [ $gwbase == snews20201116 ]; then realgweventtime="2020-11-16T14:29:01.350000Z"
  elif [ $gwbase == snews20201117 ]; then realgweventtime="2020-11-17T14:29:01.350000Z"
  elif [ $gwbase == snews20201118 ]; then realgweventtime="2020-11-18T14:29:01.350000Z"
  elif [ $gwbase == snews20201119 ]; then realgweventtime="2020-11-19T14:29:01.350000Z"
  elif [ $gwbase == snews20201120 ]; then realgweventtime="2020-11-20T14:29:01.350000Z"
  elif [ $gwbase == snews20201121 ]; then realgweventtime="2020-11-21T14:29:01.350000Z"


  else
    echo I do not know when $gwbase was.  Edit env.sh
    if ! [ "$PS1" ]; then
      exit 1
    fi
  fi

  intgwunixtime=$(TZ=UTC date +%s -d "$(echo $realgweventtime | sed -e 's/\..*//' -e 's/T/ /')")
  gwunixtime=$intgwunixtime.$(echo $realgweventtime | sed -e 's/.*\.//' -e s/Z//)
  nddef=strait-ligo-coincidence-artdaq-$gwunixtime-neardet-ligo
  fddef=strait-ligo-coincidence-artdaq-$gwunixtime-fardet-ligo

  mkdir -p $outhistdir $outhadddir

  if ! echo $PATH | grep -qE "$SRT_PRIVATE_CONTEXT/ligo([:/]|$)"; then
    PATH+=:$SRT_PRIVATE_CONTEXT/ligo
  fi
fi
