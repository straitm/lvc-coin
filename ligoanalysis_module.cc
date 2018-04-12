////////////////////////////////////////////////////////////////////////
/// \brief The ligoanalysis module looks at GW coincidences.
///
/// Given a time window, specified as an absolute time and a delta, it
/// searches for a set of interesting event types so that you can see
/// if there is a spike coincident with a gravitational wave event, or,
/// more broadly, any sort of external event you might think of. The
/// output is a set of histograms.
///
/// \author M. Strait
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"

#include "Geometry/Geometry.h"

#include "DAQDataFormats/RawTriggerTime.h"
#include "DAQDataFormats/RawEvent.h"
#include "DAQDataFormats/RawTrigger.h"
#include "DAQDataFormats/RawTriggerMask.h"
#include "DAQDataFormats/RawDataBlock.h"
#include "StandardRecord/SREnums.h"

#include "RawData/FlatDAQData.h"
#include "RawData/RawTrigger.h"

#include "RecoBase/Track.h"

#include "CelestialLocator/CelestialLocator.h"
#include "NovaTimingUtilities/TimingUtilities.h"

#include "TH2.h"

#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include <signal.h>

#include "func/timeutil.h"

// "`-._,-'"`-._,-'"`-._,-' BEGIN sky map stuff "`-._,-'"`-._,-'"`-._,-'
#include "healpix_base.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"

#include "alm.h" // Alm<T>
#include "alm_healpix_tools.h" // alm2map(), etc.
#include "alm_powspec_tools.h" // smoothWithGauss()

// The sky map from LIGO/Virgo, if available and necessary, i.e. if we
// are analyzing events with pointing
static Healpix_Map<float> * healpix_skymap = NULL;

// The critical probability density value in the sky map above which
// we are in the 90% confidence level region.  Set in beginJob().
static double skymap_crit_val = 0;
// "`-._,-'"`-._,-'"`-._,-'  END sky map stuff  "`-._,-'"`-._,-'"`-._,-'


static const int TDC_PER_US = 64;
static const int US_PER_MICROSLICE = 50; // I hope this is always true

// Set from the FCL parameters
static double gwevent_unix_double_time = 0;
static long long window_size_s = 1000;

static int gDet = caf::kUNKNOWN;

// Types of analysis, dependent on which trigger we're looking at
enum analysis_class_t { NDactivity, LiveTime, UpMu, DDenergy,
                        JustTrigger, MeV, MAX_ANALYSIS_CLASS };

class ligoanalysis : public art::EDProducer {
  public:
  explicit ligoanalysis(fhicl::ParameterSet const& pset);
  virtual ~ligoanalysis() { }; // compiles, but does not run, without this
  void produce(art::Event& evt);

  void beginJob();

  /// \brief User-supplied type of trigger being examined.
  ///
  /// Which histograms we make depends on this.
  analysis_class_t fAnalysisClass;

  /// \brief The user-supplied time of the gravitational wave burst, or
  /// whatever time you want to center the search on.
  ///
  /// Expressed in RFC-3339 format, always in UTC, always with a Z at
  /// the end (to emphasize that it is UTC).
  std::string fGWEventTime;

  /// The file name of the LIGO/Virgo skymap, or the empty string if this
  /// analysis does not care about pointing or no map is available.
  std::string fSkyMap;

  /// \brief The user-supplied length of the search window in seconds.
  ///
  /// We search for half this length on either side of the given time.
  float fWindowSize;
};

/**********************************************************************/
/*                             Histograms                             */
/**********************************************************************/

struct ligohist{
  // Things per second, for whatever it is we're selecting
  TH1D * sig = NULL;

  // Livetime per second.  This means the amount of time read out by
  // triggers in the analysis window, not the amount of time a trigger
  // was enabled in the DAQ or anything similar.
  TH1D * live = NULL;

  // Base name for the histograms.  Name for the livetime histogram will be
  // this with "live" appended.
  std::string name;

  // Answers the question "Are livetime histograms meaningful for this?"
  //
  // They aren't meaningful if the sort of event we're looking at causes
  // the trigger to fire. For instance, it doesn't mean anything to look
  // at the livetime provided by the Upmu trigger when looking at a
  // count of upward going muons. However, it may be meaningful to look
  // at the livetime of Upmu if looking at low-energy events that happen
  // to get read out because of Upmu.
  //
  // If false, we won't write out the second histogram. In some cases,
  // whether or not they are meaningful depends on the trigger that is
  // using this histogram. We'll enable this if it is at least sometimes
  // meaningful. (Or would it be better to make different histograms and
  // avoid ever having irrelevant output?)
  bool dolive;

  ligohist(const std::string & name_, const bool dolive_)
  {
    name = name_;
    dolive = dolive_;
  }
};

// Variant of ligohist with 2D histograms.
struct ligohist2d{
  TH2D * sig = NULL;

  // Binning in the second dimension
  int nbins;
  double low, high;

  // Assume livetime is the same for all detector regions, or whatever
  // we're using the second dimension for! I suppose this could be
  // violated in exceptional cases, but I hope that never matters.
  TH1D * live = NULL;

  std::string name;
  bool dolive;

  ligohist2d(const std::string & name_, const bool dolive_,
             const int nbins_, const double low_, const double high_)
  {
    name = name_;
    dolive = dolive_;
    nbins = nbins_;
    low = low_;
    high = high_;
  }
};

// Construct the histograms and hook them up with TFileService to be saved
static void init_lh(ligohist & lh)
{
  if(lh.sig != NULL || lh.live != NULL){
    fprintf(stderr, "%s already initialized.\n", lh.name.c_str());
    exit(1);
  }

  art::ServiceHandle<art::TFileService> t;

  lh.sig = t->make<TH1D>(lh.name.c_str(), "",
    window_size_s, -window_size_s/2, window_size_s/2);

  if(lh.dolive)
    lh.live = t->make<TH1D>(Form("%slive", lh.name.c_str()), "",
      window_size_s, -window_size_s/2, window_size_s/2);
}

static void init_lh2d(ligohist2d & lh)
{
  if(lh.sig != NULL || lh.live != NULL){
    fprintf(stderr, "%s already initialized.\n", lh.name.c_str());
    exit(1);
  }

  art::ServiceHandle<art::TFileService> t;

  lh.sig = t->make<TH2D>(lh.name.c_str(), "",
    window_size_s, -window_size_s/2, window_size_s/2,
    lh.nbins, lh.low, lh.high);

  if(lh.dolive)
    lh.live = t->make<TH1D>(Form("%slive", lh.name.c_str()), "",
      window_size_s, -window_size_s/2, window_size_s/2);
}

/**********************************************************************/

/*
  Get the length of the event in TDC ticks, typically 550*64, and
  "delta_tdc", the time between the trigger time and the time the event
  starts. You can subtract this off of the time that the offline gives
  each hit to get the time since the beginning of the readout, and with
  the event length, the time until the end of the readout.

  delta_tdc is a signed 64 bit integer, even though it should always be
  a small positive number, just in case. Ditto for the event length.

  Returns whether this information was successfully extracted.
*/
static bool delta_and_length(int64_t & event_length_tdc,
  int64_t & delta_tdc,
  const art::Handle< std::vector<rawdata::FlatDAQData> > & flatdaq,
  const art::Handle< std::vector<rawdata::RawTrigger> > & rawtrigger)
{
  daqdataformats::RawEvent raw;
  if(flatdaq->empty()) return false;

  raw.readData((*flatdaq)[0].getRawBufferPointer());
  if(raw.getDataBlockNumber() == 0) return false;

  raw.setFloatingDataBlock(0);
  daqdataformats::RawDataBlock& datablock = *raw.getFloatingDataBlock();

  uint64_t event_start_time = 0xffffffffffffffff;
  uint64_t event_end_time   = 0x0000000000000000;

  for(unsigned int di = 0; di < raw.getDataBlockNumber(); di++){
    raw.setFloatingDataBlock(di);
    datablock = (*raw.getFloatingDataBlock());

    if(datablock.getHeader()->getMarker() ==
         daqdataformats::datablockheader::SummaryBlock_Marker ||
       !datablock.getHeader()->checkMarker()) continue;

    for(unsigned int mi = 0; mi < datablock.getNumMicroBlocks(); mi++){
      datablock.setFloatingMicroBlock(mi);
      daqdataformats::RawMicroBlock * ub = datablock.getFloatingMicroBlock();

      // The time is always in the second and third words of the
      // microslice, which follows two words of microblock header, so
      // just get it. Justin says you can also get it from getTime(),
      // but this already works and I'm not touching it.
      const uint32_t t_marker_low  = ((uint32_t *)(ub->getBuffer()))[3];
      const uint32_t t_marker_high = ((uint32_t *)(ub->getBuffer()))[4];

      uint64_t time_marker = t_marker_low;
      time_marker |= (uint64_t)t_marker_high << 32;
      if(time_marker < event_start_time) event_start_time = time_marker;
      if(time_marker > event_end_time  ) event_end_time   = time_marker;
    }
  }

  delta_tdc = (int64_t)((*rawtrigger)[0].fTriggerTimingMarker_TimeStart
                        - event_start_time);

  // Assume that microblocks are always 50us. I hope that's true for all
  // relevant data.
  event_length_tdc = ((int64_t)(event_end_time - event_start_time))
                     + US_PER_MICROSLICE*TDC_PER_US;
  return true; // ok
}

/*********************************************************************/
/************************* Begin histograms **************************/
/*********************************************************************/

// Count of triggers, with no examination of the data within
static ligohist lh_rawtrigger("rawtrigger", false);

// Number of hits, with no filtering of any sort
static ligohist lh_rawhits("rawhits", true);

// Number of hits that are not in any Slicer4D slice, i.e. they are in the
// Slicer4D noise slice.
static ligohist lh_unslice4ddhits("unslice4ddhits", true);

// My experience with my neutron capture analysis in the ND is that
// things become well-behaved at about this level.
static const double bighit_threshold_nd = 35; // PE
static const double bighit_threshold_fd = 30;
static       double bighit_threshold    =  0; // set when we know det

// Number of unsliced hits that are over the above PE threshold
static ligohist lh_unsliced_big_hits("unslicedbighits", true);

// Number of unsliced hits that are over some PE threshold and are paired
// with a hit in an adjacent plane.  (Not an adjacent cell, because then
// we select lots of pairs from noisy modules.  We could build a noise map
// to fix that, but it would require two passes through the data.)
static ligohist lh_unsliced_hit_pairs("unslicedhitpairs", true);
static ligohist2d lh_unsliced_hit_pairs_plane
  ("unslicedhitpairsplane", true, 28, 0, 28*32);

// Count of slices with nothing around the edges, regardless of what sorts of
// objects are inside.
static ligohist lh_contained_slices("contained_slices", true);

// Number of triggers above two cuts for DDEnergy, the first pair
// in raw ADC, the second ADC per unit time.
static ligohist lh_ddenergy_locut("energy_low_cut", false);
static ligohist lh_ddenergy_hicut("energy_high_cut", false);
static ligohist lh_ddenergy_lopertime("energy_low_cut_pertime", false);
static ligohist lh_ddenergy_hipertime("energy_high_cut_pertime", false);

/********************** Histogram with pointing ***********************/

// Raw number of tracks
static ligohist lh_tracks("tracks", true);
static ligohist lh_tracks_point("tracks_point", true);

// Number of tracks with at least one contained endpoint, with some
// additional sanity checks.
static ligohist lh_halfcontained_tracks("halfcontained_tracks", true);
static ligohist lh_halfcontained_tracks_point("halfcontained_tracks_point",
                                              true);

// Number of tracks with at two contained endpoints, with some additional
// sanity checks.
static ligohist lh_fullycontained_tracks("fullycontained_tracks", true);
static ligohist lh_fullycontained_tracks_point("fullycontained_tracks_point",
                                               true);

// Number of tracks that pass the Upmu analysis
static ligohist lh_upmu_tracks("upmu_tracks", false);
static ligohist lh_upmu_tracks_point("upmu_tracks_point", false);

/*********************************************************************/
/**************************  End histograms **************************/
/*********************************************************************/

static void init_mev_hists()
{
  init_lh(lh_rawhits);
  init_lh(lh_unslice4ddhits);
  init_lh(lh_unsliced_big_hits);
  init_lh(lh_unsliced_hit_pairs);
  init_lh2d(lh_unsliced_hit_pairs_plane);
}

static void init_track_and_contained_hists()
{
  init_lh(lh_tracks);
  init_lh(lh_halfcontained_tracks);
  init_lh(lh_fullycontained_tracks);
  init_lh(lh_contained_slices);

  init_lh(lh_tracks_point);
  init_lh(lh_halfcontained_tracks_point);
  init_lh(lh_fullycontained_tracks_point);
}

// Probabilities for a point in the sky map. One is the raw value and
// the second is scaled by area, i.e. multiplied by sin(theta)
struct ac_raw{
  double raw;
  double area_corrected;
};

static bool compare_ac_raw(const ac_raw & a, const ac_raw & b)
{
  return a.raw > b.raw;
}

// Takes a skymap and smears it with a gaussian to take into account our
// detector resolution.
//
// Mostly copied from the example code in the Healpix package in
// smoothing_cxx_module.cc
//
// I'd prefer to return a smeared map from the const& argument, but I
// can't immediately see how to manage that.
static void smear_skymap(Healpix_Map<float> * map)
{
  // Some power of two between 32 and 1024, by grepping.
  // Is this a good value?
  const int nlmax = 512;

  // 3 or 5 in all examples, more often 3.
  const int num_iter = 3;

  const tsize nmod = map->replaceUndefWith0();
  if(nmod!=0)
    printf("smear_skymap() WARNING: replaced %lu undefined map pixels "
           "with a value of 0\n", nmod);

  const float avg = map->average();
  map->Add(-avg);

  // I *think* the LIGO/Virgo skymaps are not weighted, so this is
  // easier than using get_right_weights(), which requires a paramfile.
  arr<double> weight;
  weight.alloc(2*healpix_skymap->Nside());
  weight.fill(1);

  Alm<xcomplex<float> > alm(nlmax, nlmax);
  if(map->Scheme() == NEST) map->swap_scheme();

  map2alm_iter(*map, alm, num_iter, weight);

  // What really is the best resolution to use? We don't have a solid
  // number, and obviously it is a function of energy, too. This "about
  // 1 degree" is presumably 1 sigma, but it looks like smoothWithGauss
  // takes FWHM.
  const double fwhm = 1.0 /* 1 degree */ * 2.355 * M_PI/180.;

  smoothWithGauss(alm, fwhm);
  alm2map(alm, *map);

  map->Add(avg);
}

void ligoanalysis::beginJob()
{
  if(fSkyMap == "") return;

  healpix_skymap = new Healpix_Map<float>;

  try       { read_Healpix_map_from_fits(fSkyMap, *healpix_skymap); }
  catch(...){ exit(1);                                             }

  smear_skymap(healpix_skymap);

  double sumprob = 0;

  // We need to integrate the probability and find the critical value
  // above which we are in the 90% CL region. The convenient way to
  // retrieve probabilities from the map is interpolated values at
  // unevenly spaced points. Scale these by their effective area to do
  // the intergration, but set the critical value using unscaled values.
  std::vector<ac_raw> vals;

  int ni = 1000, nj = 1000;
  for(int i = 1; i < ni; i++){
    for(int j = 0; j < nj; j++){
      const double theta = i*M_PI/ni,
                   phi   = j*2.*M_PI/nj;
      const float val = healpix_skymap->interpolated_value(
        pointing(theta, // dec: except this is 0 to pi, and dec is pi/2 to -pi/2
                 phi)); // ra: as normal: 0h, 24h = 2pi
      sumprob += val * sin(theta);
      ac_raw new_ac_raw;
      new_ac_raw.raw = val;
      new_ac_raw.area_corrected = val*sin(theta);
      vals.push_back(new_ac_raw);
    }
  }

  sort(vals.begin(), vals.end(), compare_ac_raw);

  const float CL = 0.9; // 90% confidence level

  float acc = 0;
  for(unsigned int i = 0; i < vals.size(); i++){
    acc += vals[i].area_corrected/sumprob;
    if(acc > CL){
      skymap_crit_val = vals[i].raw;
      break;
    }
  }

  // Print the map to the screen just so we know something is happening
  printf("Sky map %.0f%% region:\n", CL*100);
  const int across = 80, down = 40;
  for(int i = 0; i < down; i++){
    for(int j = across-1; j >= 0; j--){
      const double theta = i*M_PI/down, phi = j*2*M_PI/across;
      const float val = healpix_skymap->interpolated_value(pointing(theta, phi));
      printf("%c", val > skymap_crit_val?'X':'-');
    }
    printf("\n");
  }
}

ligoanalysis::ligoanalysis(fhicl::ParameterSet const& pset) : EDProducer(),
  fGWEventTime(pset.get<std::string>("GWEventTime")),
  fSkyMap(pset.get<std::string>("SkyMap")),
  fWindowSize(pset.get<unsigned long long>("WindowSize"))
{
  const std::string analysis_class_string(pset.get<std::string>("AnalysisClass"));

  if     (analysis_class_string == "NDactivity") fAnalysisClass = NDactivity;
  else if(analysis_class_string == "LiveTime")   fAnalysisClass = LiveTime;
  else if(analysis_class_string == "UpMu")       fAnalysisClass = UpMu;
  else if(analysis_class_string == "DDenergy")   fAnalysisClass = DDenergy;
  else if(analysis_class_string == "JustTrigger")fAnalysisClass = JustTrigger;
  else if(analysis_class_string == "MeV")        fAnalysisClass = MeV;
  else{
    fprintf(stderr, "Unknown AnalysisClass \"%s\" in job fcl. See list "
            "in ligoanalysis.fcl.\n", analysis_class_string.c_str());
    exit(1);
  }

  gwevent_unix_double_time = rfc3339_to_unix_double(fGWEventTime);
  window_size_s = fWindowSize;

  switch(fAnalysisClass){
    case NDactivity:
      init_lh(lh_rawtrigger);
      init_track_and_contained_hists();
      break;
    case JustTrigger:
      init_lh(lh_rawtrigger);
      break;
    case UpMu:
      init_lh(lh_upmu_tracks);
      init_lh(lh_upmu_tracks_point);
      break;
    case LiveTime:
      init_mev_hists();
      init_track_and_contained_hists();
      break;
    case DDenergy:
      init_lh(lh_rawtrigger);
      init_lh(lh_ddenergy_locut);
      init_lh(lh_ddenergy_hicut);
      init_lh(lh_ddenergy_lopertime);
      init_lh(lh_ddenergy_hipertime);
      break;
    case MeV:
      init_mev_hists();
      break;
    default:
      printf("No case for type %d\n", fAnalysisClass);
  }
}

/**********************************************************************/
/*                          The meat follows                          */
/**********************************************************************/

// Get the FlatDAQData, either from "minbias", in the case of supernova MC overlays,
// or "daq", for everything else.
static void getflatdaq(art::Handle< std::vector<rawdata::FlatDAQData> > & flatdaq,
                       const art::Event & evt)
{
  evt.getByLabel("minbias", flatdaq);
  if(flatdaq.failedToGet())
    evt.getByLabel("daq", flatdaq);
}

static void getrawtrigger(art::Handle< std::vector<rawdata::RawTrigger> > & trg,
                         const art::Event & evt)
{
  evt.getByLabel("minbias", trg);
  if(trg.failedToGet())
    evt.getByLabel("daq", trg);
}

static void getrawdigits(art::Handle< std::vector<rawdata::RawDigit> > & digits,
                         const art::Event & evt)
{
  evt.getByLabel("minbias", digits);
  if(digits.failedToGet())
    evt.getByLabel("daq", digits);
}

// Return the livetime in this event in seconds, as it is relevant for
// raw hits (same as for anything else? Maybe not if we don't trust
// tracks close to the time-edges of events).
static double rawlivetime(const art::Event & evt)
{
  art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
  getflatdaq(flatdaq, evt);
  if(flatdaq.failedToGet()){
    static bool first = true;
    if(first)
      puts("No FlatDAQ. Probably unoverlayed MC. Returning zero live time");
    first = false;
    return 0;
  }

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  getrawtrigger(rawtrigger, evt);

  if(rawtrigger->empty()) return -1;

  int64_t event_length_tdc = 0, delta_tdc = 0;

  if(!delta_and_length(event_length_tdc, delta_tdc, flatdaq, rawtrigger))
    return -1;

  return event_length_tdc / TDC_PER_US * 1e-6;
}

// Add 'sig' to the signal and 'live' to the livetime in bin number 'bin'.
static void THplusequals(ligohist & lh, const int bin, const double sig,
                         const double live)
{
  // Use SetBinContent instead of Fill(x, weight) to avoid having to look up
  // the bin number twice.
  lh.sig->SetBinContent(bin, lh.sig->GetBinContent(bin) + sig);
  if(lh.dolive)
    lh.live->SetBinContent(bin, lh.live->GetBinContent(bin) + live);
}

// Same as THplusequals for a ligohist2d. Note inconsistent arguments.
// It wants the bin number for time, but the value of the variable for
// the second dimension.
static void THplusequals2d(ligohist2d & lh, const int timebin,
                           const double othervalue, const double sig,
                           const double live)
{
  const int otherbin = lh.sig->GetYaxis()->FindBin(othervalue);
  lh.sig->SetBinContent(timebin, otherbin,
                        lh.sig->GetBinContent(timebin, otherbin) + sig);
  if(lh.dolive)
    lh.live->SetBinContent(timebin, lh.live->GetBinContent(timebin) + live);
}

// Return the bin number for this event, i.e. the number of seconds from the
// beginning of the window, plus 1.
static int timebin(const art::Event & evt)
{
  const double evt_time = art_time_to_unix_double(evt.time().value());
  return floor(evt_time - gwevent_unix_double_time) + window_size_s/2
         + 1; // stupid ROOT 1-based numbering!
}

// Returns true if a point is "contained" for purposes of deciding if a track
// or physics slice is contained.
static bool contained(const TVector3 & v)
{
  if(gDet == caf::kNEARDET)
    return fabs(v.X()) < 150 && fabs(v.Y()) < 150 && v.Z() > 40 && v.Z() < 1225;
  if(gDet == caf::kFARDET)
    return fabs(v.X()) < 650 &&
      v.Y() < 500 && v.Y() > -650 &&
      v.Z() > 75 && v.Z() < 5900;
  fprintf(stderr, "Unknown detector %d\n", gDet);
  exit(1);
}

// Returns true if the track enters and exits.  Or if both of its ends
// are pretty close to the edge.  But not if it is just steep.
static bool un_contained_track(const rb::Track & t)
{
  return !contained(t.Start()) && !contained(t.Stop());
}

// Check that the reconstruction (probably BreakPointFitter) doesn't think
// either end of the track points nearly along a plane AND check that if we
// just draw a line from one end of the track to the other, that that doesn't
// either.  This second part is to make sure that some tiny kink at the start
// or end can't fool us into thinking it is well-contained.
static bool good_track_direction(const rb::Track & t)
{
  static const double tan_track_cut = atan(15 * M_PI/180);

  const double tot_dx = t.Start().X() - t.Stop().X();
  const double tot_dy = t.Start().Y() - t.Stop().Y();
  const double tot_dz = t.Start().Z() - t.Stop().Z();

  const double rec_dx1 = t.    Dir().X();
  const double rec_dx2 = t.StopDir().X();
  const double rec_dy1 = t.    Dir().Y();
  const double rec_dy2 = t.StopDir().Y();
  const double rec_dz1 = t.    Dir().Z();
  const double rec_dz2 = t.StopDir().Z();

  return fabs(tot_dz /tot_dx ) > tan_track_cut &&
         fabs(tot_dz /tot_dy ) > tan_track_cut &&
         fabs(rec_dz1/rec_dx1) > tan_track_cut &&
         fabs(rec_dz2/rec_dx2) > tan_track_cut &&
         fabs(rec_dz1/rec_dy1) > tan_track_cut &&
         fabs(rec_dz2/rec_dy2) > tan_track_cut;
}

// Fewest planes we'll accept a track having, or in the case of fully-contained
// tracks, the fewest *contiguous* planes.
static const int min_plane_extent = 10;

// Returns true if the track starts AND stops inside the detector,
// and we're pretty sure that this isn't artefactual.
static bool fully_contained_track(const rb::Track & t)
{
  // Note how here we're looking at all the hits' positions, not
  // just the line whose endpoints are Start() and Stop().  For steep
  // tracks, this reveals their uncontainedness.
  return contained(t.MaxXYZ()) && contained(t.MinXYZ()) &&
          good_track_direction(t) &&

          // Don't accept a track just because of a stray hit that gets
          // swept into it.
          t.MostContiguousPlanes(geo::kXorY) >= min_plane_extent;
}

// Returns true if either end of the track is uncontained or might be
static bool half_uncontained_track(const rb::Track & t)
{
  return !contained(t.Start()) || !contained(t.Stop()) ||
          !good_track_direction(t) ||
          t.ExtentPlane() < min_plane_extent;
}

// Returns true if the track either starts or stops in the detector, or both,
// and we're pretty sure that this isn't artefactual.
static bool half_contained_track(const rb::Track & t)
{
  return (contained(t.Start()) || contained(t.Stop())) &&
          good_track_direction(t) &&
          t.ExtentPlane() >= min_plane_extent;
}

// return the index in the slice array that the given track is in
static int which_slice_is_this_track_in(
  const rb::Track & t,
  const art::Handle< std::vector<rb::Cluster> > & slice)
{
  // I'm sure this is not the best way, but I have had it with trying to
  // figure out what the best way is, and I'm just going to do it *some*
  // way.
  if(t.NCell() == 0) return -1;
  const art::Ptr<rb::CellHit> ahit =
    t.Cell(0); // some random hit on the track

  // Skip slice 0 since it is the noise slice
  for(unsigned int i = 1; i < slice->size(); i++){
    const rb::Cluster & slc = (*slice)[i];
    for(unsigned int j = 0; j < slc.NCell(); j++){
      const art::Ptr<rb::CellHit> shit = slc.Cell(j);
      if(*ahit == *shit) return i;
    }
  }
  return -1;
}

/************************************************/
/* Here come the functions that fill histograms */
/************************************************/

// Put put raw triggers into a histogram
static void count_triggers(const art::Event & evt)
{
  THplusequals(lh_rawtrigger, timebin(evt), 1, 1);
}

static void count_ddenergy(const art::Event & evt)
{
  art::Handle< std::vector<rawdata::RawDigit> > rawhits;
  getrawdigits(rawhits, evt);

  int64_t sumadc = 0;

  for(unsigned int i = 0; i < rawhits->size(); i++)
    sumadc += (*rawhits)[i].ADC();

  const double rawtime = rawlivetime(evt);

  printf("ADC: %ld %12.0f\n", sumadc, sumadc/rawtime);

  if(sumadc/rawtime > 5e10)
    THplusequals(lh_ddenergy_lopertime, timebin(evt), 1, rawtime);
  if(sumadc/rawtime > 5e11)
    THplusequals(lh_ddenergy_hipertime, timebin(evt), 1, rawtime);

  if(sumadc >  5000000)
    THplusequals(lh_ddenergy_locut, timebin(evt), 1, rawtime);
  if(sumadc > 50000000)
    THplusequals(lh_ddenergy_hicut, timebin(evt), 1, rawtime);
}

// Count the number of raw hits in the event and fill the appropriate histograms
static void count_hits(const art::Event & evt)
{
  art::Handle< std::vector<rawdata::RawDigit> > rawhits;
  getrawdigits(rawhits, evt);

  printf("Hits in this event: %lu\n", rawhits->size());

  THplusequals(lh_rawhits, timebin(evt), rawhits->size(), rawlivetime(evt));
}

// Minimal hit information, distilled from CellHit
struct mhit{
  float tns; // fine timing, in nanoseconds
  float tpos; // transverse position
  uint16_t plane;
  int cell;
  bool used; // has this hit been used in a pair yet?
};

// Minimal slice information, distilled from rb::Cluster
struct mslice{
  float mintns, maxtns;
  uint16_t minplane, maxplane;
  uint16_t mincellx, maxcellx, mincelly, maxcelly;
};

// Helper function for count_unsliced_hit_pairs().
// Builds list of distilled slice information
static std::vector<mslice> make_sliceinfo_list(
  const art::Handle< std::vector<rb::Cluster> > & slice)
{
  std::vector<mslice> sliceinfo;

  // Start with slice number 1 because 0 is the "noise slice".
  for(unsigned int i = 1; i < slice->size(); i++){
    mslice slc;
    slc.mintns = (*slice)[i].MinTNS();
    slc.maxtns = (*slice)[i].MaxTNS();
    slc.minplane = (*slice)[i].MinPlane();
    slc.maxplane = (*slice)[i].MaxPlane();
    slc.mincellx = (*slice)[i].MinCell(geo::kX);
    slc.maxcellx = (*slice)[i].MaxCell(geo::kX);
    slc.mincelly = (*slice)[i].MinCell(geo::kY);
    slc.maxcelly = (*slice)[i].MaxCell(geo::kY);
    sliceinfo.push_back(slc);
  }
  return sliceinfo;
}

// Helper function for count_unsliced_hit_pairs().  Selects hits that
// are candidates to be put into hit pairs.
static std::vector<mhit> select_hits_for_mev_search(
  const rb::Cluster & noiseslice, const std::vector<mslice> & sliceinfo)
{
  art::ServiceHandle<geo::Geometry> geo;

  std::vector<mhit> mhits;


  for(unsigned int i = 0; i < noiseslice.NCell(); i++){
    const float fd_low_adc =  85, fd_high_adc =  600;
    const float nd_low_adc = 107, nd_high_adc = 2500;

    const float low_adc  = gDet == caf::kNEARDET? nd_low_adc : fd_low_adc;
    const float high_adc = gDet == caf::kNEARDET? nd_high_adc: fd_high_adc;

    if(noiseslice.Cell(i)->ADC() <=  low_adc) continue;
    if(noiseslice.Cell(i)->ADC() >= high_adc) continue;

    const int cell  = noiseslice.Cell(i)->Cell();
    const int plane = noiseslice.Cell(i)->Plane();

    if(gDet == caf::kNEARDET){
      const int ndnplaneedge = 4;

      // Exclude the whole muon catcher. Can't reasonably have a
      // supernova-type event hit planes on each side of a steel plane,
      // and while they might hit adjacent scintillator planes in the
      // muon catcher, I don't want to deal with all the additional
      // complications there.
      if(plane <= ndnplaneedge) continue;
      if(plane >= 192-ndnplaneedge) continue;

      // At top of detector, stricter cut based on observed backgrounds
      if(noiseslice.Cell(i)->View() == geo::kY && cell >= 96 - 20) continue;

      // Drop outermost cells. This excludes most muon tracks that
      // just barely enter the detector, but don't get reconstructed.
      const int ndncelledge_x = 4;
      if(noiseslice.Cell(i)->View() == geo::kX &&
         (cell <= ndncelledge_x || cell >= 96 - ndncelledge_x)) continue;
    }
    else{ // far
      const int fdnplaneedge = 2;
      if(plane <= fdnplaneedge) continue;
      if(plane >= 895-fdnplaneedge) continue;

      if(noiseslice.Cell(i)->View() == geo::kY && cell >= 383 - 50)
        continue;

      const int fdncelledge_x = 10;
      if(noiseslice.Cell(i)->View() == geo::kX &&
         (cell <= fdncelledge_x || cell >= 383 - fdncelledge_x)) continue;
    }

    const float tns = noiseslice.Cell(i)->TNS();

    bool pass = true;

    // Exclude a rectangular box around each slice for a given time
    // period before and after the slice, with a given spatial buffer.
    for(unsigned int j = 0; j < sliceinfo.size(); j++){
      // Optimized for the FD.  The ND isn't sensitive to these.
      const float time_until_slc_cut =  2e3;
      const float time_since_slc_cut = 13e3;

      // Geometrically about correct, but perhaps should be scaled by density
      // or radiation length or neutron cross section or something. Or not,
      // since which of those is right depends on what you're looking at.
      const double planes_per_cell = 76./39.;

      // Optimized for the FD, and at the ND, you can set this to any
      // positive value and it has the same result.
      const int planebuffer = 15;
      const int cellbuffer = planebuffer * planes_per_cell;

      if(tns > sliceinfo[j].mintns - time_until_slc_cut &&
         tns < sliceinfo[j].maxtns + time_since_slc_cut &&
         plane > sliceinfo[j].minplane - planebuffer &&
         plane < sliceinfo[j].maxplane + planebuffer){

        if(noiseslice.Cell(i)->View() == geo::kX){
          if(cell > sliceinfo[j].mincellx - cellbuffer &&
             cell < sliceinfo[j].maxcellx + cellbuffer) pass = false;
        }
        else{ // Y view
          if(cell > sliceinfo[j].mincelly - cellbuffer &&
             cell < sliceinfo[j].maxcelly + cellbuffer) pass = false;
        }
      }
    }
    if(!pass) continue;

    mhit h;
    h.used  = false;
    h.tns   = tns;
    h.plane = plane;
    h.cell  = cell;
    h.tpos  = geo->CellTpos(plane, cell);

    mhits.push_back(h);
  }
  return mhits;
}

static void count_unsliced_hit_pairs(const art::Event & evt)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);
  if(slice.failedToGet()){
    printf("Unexpected lack of slicer product!\n");
    return;
  }

  if(slice->empty()){
    printf("Unexpected event with zero slices!\n");
    return;
  }

  // Make list of all times and locations of "physics" slices. Exclude
  // other hits near them in time to drop Michels & maybe also neutrons.
  const std::vector<mslice> sliceinfo = make_sliceinfo_list(slice);

  // Find hits which we'll accept for possible membership in pairs.
  // Return value is not const because we modify ::used below.
  std::vector<mhit> mhits = select_hits_for_mev_search((*slice)[0], sliceinfo);

  // Intentionally bigger than optimum value suggested by MC (150ns FD,
  // 10ns(!) ND) because we know that the data has a bigger time spread
  // from, for instance, the very well known slice duration discrepancy
  // (see, e.g., doc-19053).
  const float timewindow = 250; // ns

  unsigned int hitpairs = 0;

  // Rough effective light speed in fiber
  const float lightspeed = 17.; // cm/ns

  for(unsigned int i = 0; i < mhits.size(); i++){
    for(unsigned int j = 0; j < mhits.size(); j++){
      if(i == j) continue;
      if(mhits[i].used || mhits[j].used) continue;

      if(abs(mhits[j].plane - mhits[i].plane) != 1) continue;

      const float time_i_corr = mhits[i].tns + mhits[j].tpos / lightspeed;
      const float time_j_corr = mhits[j].tns + mhits[i].tpos / lightspeed;

      if(fabs(time_j_corr - time_i_corr) > timewindow) continue;

      hitpairs++;

      THplusequals2d(lh_unsliced_hit_pairs_plane, timebin(evt), mhits[i].plane,
                     1, rawlivetime(evt));

      #ifdef PRINTCELL
        const int eveni = j%2 == 0?j:i,
                  oddi  = j%2 == 1?j:i;
        printf("Bighitpair: %4d %4d  %4d %4d\n",
               mhits[eveni].plane, mhits[eveni].cell, // y
               mhits[oddi ].plane, mhits[oddi ].cell); // x
      #endif

      // Do not allow these hits to be used again. This makes the number
      // of pairs found slightly dependent on the order of the original
      // list, for instance if the middle two of a group of four are
      // picked first, you get one pair instead of two, but I think it
      // is not a big deal.
      mhits[i].used = mhits[j].used = true;
      break;
    }
  }

  printf("Big hit pairs in the noise slice: %u\n", hitpairs);

  THplusequals(lh_unsliced_hit_pairs, timebin(evt), hitpairs,
               rawlivetime(evt));
}

static void count_unsliced_big_hits(const art::Event & evt)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);
  if(slice.failedToGet()){
    printf("Unexpected lack of slicer product\n");
    return;
  }

  if(slice->empty()){
    printf("Unexpected event with zero slices!\n");
    return;
  }

  unsigned int bighits = 0;

//#define PRINTPE

  for(unsigned int i = 0; i < (*slice)[0].NCell(); i++){
    bighits += (*slice)[0].Cell(i)->PE() > bighit_threshold;
#ifdef PRINTPE
    printf("unslicedPE: %f\n", (*slice)[0].Cell(i)->PE());
#endif
  }

#ifdef PRINTPE
  for(unsigned int j = 1; j < slice->size(); j++)
    for(unsigned int i = 0; i < (*slice)[i].NCell(); i++)
      printf("slicedPE: %f\n", (*slice)[i].Cell(i)->PE());
#endif

  printf("Big hits in this event in the noise slice: %u\n", bighits);

  THplusequals(lh_unsliced_big_hits, timebin(evt), bighits,
               rawlivetime(evt));
}

// Counts hits in the noise slices and adds the results to the output histograms
static void count_unslice4dd_hits(const art::Event & evt)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);
  if(slice.failedToGet()){
    printf("Unexpected lack of slicer product!\n");
    return;
  }

  if(slice->empty()){
    printf("Unexpected event with zero slices!\n");
    return;
  }

  printf("Hits in this event in the Slicer4D noise slice: %d\n",
         (*slice)[0].NCell());

  THplusequals(lh_unslice4ddhits, timebin(evt), (*slice)[0].NCell(),
              rawlivetime(evt));
}

static void count_upmu(const art::Event & evt)
{
  art::Handle< std::vector<rb::Track> > upmu;
  evt.getByLabel("upmuanalysis", upmu);
  if(upmu.failedToGet()){
    printf("No UpMu product to read\n");
    _exit(1);
  }

  THplusequals(lh_upmu_tracks, timebin(evt),
               upmu->size(), rawlivetime(evt));
}

// Counts tracks with various cuts and contained slices and adds the
// results to the output histograms
static void count_tracks_containedslices(const art::Event & evt)
{
  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);
  if(slice.failedToGet()){
    printf("Unexpected lack of slicer product!\n");
    return;
  }

  art::Handle< std::vector<rb::Track> > tracks;
  evt.getByLabel("breakpoint", tracks);
  if(tracks.failedToGet()){
    printf("Unexpected lack of breakpoint product!\n");
    return;
  }

  const int timbin = timebin(evt);
  const double livetime = rawlivetime(evt);

  printf("Tracks in this event: %lu\n", tracks->size());
  THplusequals(lh_tracks, timbin, tracks->size(), livetime);

  // Count tracks with either a contained start or end, but not if they are
  // very steep in x or y, since those probably aren't really contained and may
  // not even be complete. Also don't count more than one per slice, nor count
  // slices with uncontained tracks.  This protects against counting a brem as
  // a contained track.

  // Find out what slice each track is in and make a list of slices with tracks
  // that aren't fully contained
  std::vector<int> trk_slices(tracks->size(), -1);
  std::set<int> slices_with_uc_tracks, slices_with_huc_tracks;
  std::set<int> contained_slices;
  std::set<int> contained_shower_slices;
  for(unsigned int i = 0; i < slice->size(); i++){
    if(!(*slice)[i].NCell(geo::kX) || !(*slice)[i].NCell(geo::kY))
      continue;

    if(!contained((*slice)[i].MinXYZ()) ||
       !contained((*slice)[i].MaxXYZ())) continue;

    contained_slices.insert(i);

    if((*slice)[i].NCell() < 10) continue;

    const int planesx = ((*slice)[i].ExtentPlane(geo::kX)+1)/2;
    const int planesy = ((*slice)[i].ExtentPlane(geo::kY)+1)/2;

    if(planesx < 5 || planesy < 5) continue;

    const int cellsx = (*slice)[i].ExtentCell(geo::kX);
    const int cellsy = (*slice)[i].ExtentCell(geo::kY);

    // Must be somewhat extended in z (not a straight-down track), and
    // whatever box is defined by the cell and plane extent must neither
    // be too sparse (random hits plus a straight-down track) nor too full
    // (FEB flashing).  A block of flashers plus a random hit is not
    // really excludable with this cut.
    //
    // These numbers are reasonable, but not very rigorous.
    const double min_boxupancy = 0.02,
                 max_boxupancy = 0.10;
    const double boxupancy_x = (*slice)[i].NXCell()/double(cellsx*planesx),
                 boxupancy_y = (*slice)[i].NYCell()/double(cellsy*planesy);

    if(boxupancy_x < max_boxupancy && boxupancy_y < max_boxupancy &&
       boxupancy_x > min_boxupancy && boxupancy_y > min_boxupancy)
      contained_shower_slices.insert(i);
  }

  // After next loop, these are true if the corresponding track in 'tracks'
  // points back to the sky region from LIGO/Virgo
  std::vector<bool> track_point;

  for(unsigned int i = 0; i < tracks->size(); i++){
    trk_slices[i] = which_slice_is_this_track_in((*tracks)[i], slice);
    if(un_contained_track((*tracks)[i]))
      slices_with_uc_tracks.insert(trk_slices[i]);
    if(half_uncontained_track((*tracks)[i]))
      slices_with_huc_tracks.insert(trk_slices[i]);

    bool track_points_to_event = false;

    if(healpix_skymap != NULL){
      double ra, dec;
      art::ServiceHandle<locator::CelestialLocator> fSunPos;

      art::Handle<std::vector<rawdata::RawTrigger> > rawtrigger;
      getrawtrigger(rawtrigger, evt);
      if(rawtrigger->empty()) return;
      unsigned long long event_time =
        (*rawtrigger)[0].fTriggerTimingMarker_TimeStart;
      struct timespec ts;
      novadaq::timeutils::
        convertNovaTimeToUnixTime(event_time, ts);

      fSunPos->GetTrackRaDec((*tracks)[i].Dir(), ts.tv_sec, ra, dec);

      // XXX double-check two things: Does pointing accept angles
      // out of range and wrap around as expected?  And is its definition
      // of declination, a.k.a. theta, really off by pi/2?
      // Accept tracks that point in the reverse direction, since the
      // orientation of the track is arbitrary.
      track_points_to_event = skymap_crit_val < std::max(
        healpix_skymap->interpolated_value(pointing( dec + M_PI_2, ra       )),
        healpix_skymap->interpolated_value(pointing(-dec + M_PI_2, ra + M_PI)));
    }
    track_point.push_back(track_points_to_event);
  }

  unsigned int ntracks_that_point = 0;
  for(unsigned int i = 0; i < tracks->size(); i++)
    if(track_point[i]) ntracks_that_point++;
  THplusequals(lh_tracks_point, timbin, ntracks_that_point, livetime);

  // Find any tracks that are half-contained/fully-contained and do not share a
  // slice with any tracks that are, like, totally uncontained, man
  std::set<int> slices_with_hc_tracks, slices_with_fc_tracks;

  // Same, but the tracks must point towards the LIGO/Virgo event
  std::set<int> slices_with_hc_tracks_point, slices_with_fc_tracks_point;

  for(unsigned int i = 0; i < tracks->size(); i++){
    // Exclude 2D tracks
    if((*tracks)[i].Stop().X() == 0 || (*tracks)[i].Stop().Y() == 0) continue;

    // To be called a slice with half-contained tracks, it must not have any
    // tracks that both enter and exit
    if(!slices_with_uc_tracks.count(trk_slices[i]) &&
       half_contained_track((*tracks)[i])){
      slices_with_hc_tracks.insert(trk_slices[i]);
      if(track_point[i])
        slices_with_hc_tracks_point.insert(trk_slices[i]);
    }

    // To be called a slice with fully contained tracks, it must not have any
    // track that either enters or exits, and the slice itself must be
    // contained.  So, in other words, no hits around the edges, and no tracks
    // reconstructed to be near the edges even in the absence of hits.
    if(!slices_with_huc_tracks.count(trk_slices[i]) &&
       fully_contained_track((*tracks)[i]) &&
       contained_slices.count(trk_slices[i])){
      slices_with_fc_tracks.insert(trk_slices[i]);
      if(track_point[i])
        slices_with_fc_tracks_point.insert(trk_slices[i]);
    }
  }

  printf("Slices with half-contained tracks: %2lu (%lu pointing)\n",
         slices_with_hc_tracks.size(), slices_with_hc_tracks_point.size());

  printf("Slices with fully-contained tracks: %2lu (%lu pointing)\n",
         slices_with_fc_tracks.size(), slices_with_fc_tracks_point.size());
  for(std::set<int>::iterator i = slices_with_fc_tracks.begin();
      i != slices_with_fc_tracks.end(); i++)
    printf("  slice %3d\n", *i);

  printf("Contained GeV physics-like slices: %lu\n",
         contained_shower_slices.size());
  for(std::set<int>::iterator i = contained_shower_slices.begin();
      i != contained_shower_slices.end(); i++)
    printf("  slice %3d\n", *i);

  THplusequals(
    lh_halfcontained_tracks, timbin, slices_with_hc_tracks.size(), livetime);
  THplusequals(
    lh_fullycontained_tracks, timbin, slices_with_fc_tracks.size(), livetime);
  THplusequals(
    lh_halfcontained_tracks_point, timbin, slices_with_hc_tracks_point.size(),
    livetime);
  THplusequals(
    lh_fullycontained_tracks_point, timbin, slices_with_fc_tracks_point.size(),
    livetime);
  THplusequals(
    lh_contained_slices, timbin, contained_shower_slices.size(), livetime);
}

static void count_mev(art::Event & evt)
{
  count_hits(evt);
  count_unslice4dd_hits(evt);
  count_unsliced_big_hits(evt);
  count_unsliced_hit_pairs(evt);
}

// TODO Somehow deal with overlapping triggers?  I found a case where a NuMI
// spill caused three overlapping ddactivity1 triggers, each offset slightly.
// And I noticed it because RemoveBeamSpills fails to remove it.  Ideally, I'd
// drop all the hits that were seen before and only reconstruct the new ones,
// but that is non-trivial.
void ligoanalysis::produce(art::Event & evt)
{
  // Must be called on every event to prevent SIGPIPE loops.
  signal(SIGPIPE, SIG_DFL);

  art::ServiceHandle<geo::Geometry> geo;
  gDet = geo->DetId();

  bighit_threshold = gDet == caf::kNEARDET?
    bighit_threshold_nd : bighit_threshold_fd;

  switch(fAnalysisClass){
    case NDactivity:
      count_triggers(evt);
      count_tracks_containedslices(evt);
      break;
    case JustTrigger:
      count_triggers(evt);
      break;
    case UpMu:
      count_upmu(evt);
      break;
    case LiveTime:
      count_mev(evt);
      count_tracks_containedslices(evt);
      break;
    case DDenergy:
      count_triggers(evt);
      count_ddenergy(evt);
      break;
    case MeV:
      count_mev(evt);
      break;
    default:
      printf("No case for type %d\n", fAnalysisClass);
  }
}

DEFINE_ART_MODULE(ligoanalysis)
