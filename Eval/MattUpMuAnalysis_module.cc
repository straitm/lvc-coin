////////////////////////////////////////////////////////////////////////
// Class:       MattUpMuAnalysis
// Module Type: producer
// File:        MattUpMuAnalysis_module.cc
//
// Generated at Mon May 16 20:16:31 2016 by Aristeidis Tsaris using artmod
// from cetpkgsupport v1_08_07.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

#include "Utilities/func/MathUtil.h"
#include "RecoBase/Track.h"
#include "RecoBase/RecoHit.h"
#include "CelestialLocator/CelestialLocator.h"
#include "RawData/DAQHeader.h"
#include "RawData/RawTrigger.h"
#include "NovaTimingUtilities/TimingUtilities.h"
#include "TrackFit/CosmicTrackUtilities.h"

#include <boost/math/special_functions.hpp> // for llr

#include "TVector3.h"
#include "TNtuple.h"
#include "TH1.h"

namespace upmuana {
  class MattUpMuAnalysis;

  inline double getErr(double PE) {
    return 165143/(1882.9+pow(PE,2.11447)) + 10.4321;
  }

  inline double getDist(double x1, double y1, double z1,
                        double x2, double y2, double z2) {
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  }

  bool recoHitComp (rb::RecoHit lhs, rb::RecoHit rhs) {
    if (lhs.Y() == rhs.Y()) {
      if (lhs.Z() == rhs.Z()) {
        return lhs.T() < rhs.T();
      }
      else return lhs.Z() < rhs.Z();
    }
    else return lhs.Y() < rhs.Y();
  }

  bool pairHitComp (std::pair<rb::CellHit, rb::RecoHit> lhp,
                    std::pair<rb::CellHit, rb::RecoHit> rhp) {
    return recoHitComp(lhp.second, rhp.second);
  }
}

class upmuana::MattUpMuAnalysis : public art::EDProducer {
public:
  explicit MattUpMuAnalysis(fhicl::ParameterSet const & p);
  void produce(art::Event & e) override;
  void beginJob() override;

private:

  art::InputTag trackTag_;
  art::InputTag fSliceModuleLabel;
  art::InputTag fSliceInstance;

  float containmentType(rb::Track);
  double fMinY = -700.0, fMinX = -650.0, fMinZ = 100.0,
    fMaxY = 400.0, fMaxX = 650.0, fMaxZ = 5800.0;

  double getLLR(std::set< std::pair<rb::CellHit,rb::RecoHit>,
		bool(*)(std::pair<rb::CellHit,rb::RecoHit>,
			std::pair<rb::CellHit,rb::RecoHit>)> ,
                std::vector<rb::RecoHit> &outliers,
                double &P_up, double &P_dn, double &chi2, double &slope);

  void LLR(std::vector<double>& eT,
           std::vector<double>& mT,
           std::vector<double>& mTErr, double& slope, double& chi2,
           double& P_up, double& P_dn, std::vector<rb::RecoHit> &in,
           std::vector<rb::RecoHit> &outliers);

  void LinFit(const std::vector<double>& x,
              const std::vector<double>& y, double *fitpar);
  void LinFit(const std::vector<double>& x,
              const std::vector<double>& y,
              const std::vector<double>& ye, double *fitpar);

  TVector3 AnglesToVector(double zen, double azi) const;

  art::ServiceHandle<locator::CelestialLocator> fSunPos;
  trk::CosmicTrackUtilities   fUtilities;

  TNtuple* ntp_track;
};

upmuana::MattUpMuAnalysis::MattUpMuAnalysis(fhicl::ParameterSet const & p)
  :
  EDProducer()
  , trackTag_(p.get<std::string>("trackInputTag"))
  , fSliceModuleLabel    (p.get<std::string>("SliceModuleLabel"))
  , fSliceInstance       (p.get<std::string>("SliceInstanceLabel"))
  , ntp_track(nullptr)
{
  produces< std::vector<rb::Track> >();
}

void upmuana::MattUpMuAnalysis::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  ntp_track = tfs->make<TNtuple>( "ntp_track", "Track Ntuple", "Run:SubRun:Event:SliceID:TrackID:Nhits:NRecohits:NOutliers:ProbUp:ProbDn:LLR:Chi2:Slope:LLRX:Chi2X:SlopeX:LLRY:Chi2Y:SlopeY:R2X:R2Y:StartX:StartY:StartZ:StartT:EndX:EndY:EndZ:EndT:TrackHitsX:TrackHitsY:Length:dirX:dirY:dirZ:eleAngle:totalE:containment:avgTX:avgTY:totalMSlices:SunZen:SunAzi:CosTheta:Azimuth:dotSun:zen_trk:azi_trk:event_time:chi2X_fit:chi2Y_fit:trackX_fromfit:trackY_fromfit:trackX_fromfit_slice:trackY_fromfit_slice");
}

void upmuana::MattUpMuAnalysis::produce(art::Event & e)
{
  // Get all slices
  art::Handle< std::vector<rb::Cluster > > slicecol;
  e.getByLabel( fSliceModuleLabel, slicecol);
  art::PtrVector<rb::Cluster > slces;

  for(int i = 0; i < (int)slicecol->size(); i++){
    art::Ptr<rb::Cluster> sli(slicecol, i);
    slces.push_back(sli);
  }

  art::Handle<std::vector<rb::Track> > tracks;
  e.getByLabel(trackTag_, tracks);
  art::FindOneP<rb::Cluster> fmSlice(tracks, e, trackTag_);

  art::Handle<std::vector<rawdata::RawTrigger> > raw_triggers;
  e.getByLabel("daq", raw_triggers); // this is need it
  if (raw_triggers->size() != 1)
    throw art::Exception(art::errors::DataCorruption)
      << __FILE__ << ":" << __LINE__ << "\n"
      << "RawTrigger vector has incorrect size: " << raw_triggers->size()
      << std::endl;

  std::vector<const rb::Cluster*> slices;
  std::vector<float> slice_containment;

  art::Handle<rawdata::DAQHeader> raw_daqheader;
  e.getByLabel("daq", raw_daqheader);

  double zen, azi;
  unsigned long long event_time(raw_triggers->front().fTriggerTimingMarker_TimeStart);

  struct timespec ts;
  novadaq::timeutils::
    convertNovaTimeToUnixTime(event_time, ts);

  fSunPos->GetSunPosition_FD(ts.tv_sec, zen, azi);

  // Save passed track to the output product
  std::unique_ptr< std::vector<rb::Track> > output(new std::vector<rb::Track>);

  std::vector<unsigned> tracks_per_slice;
  for (size_t i_track=0; i_track < tracks->size(); ++i_track){

    rb::Track theTrack = tracks->at(i_track);

    double ra, dec;
    fSunPos->GetTrackRaDec(theTrack.Dir(), ts.tv_sec, ra, dec);

    float containment = containmentType(theTrack);

    float slice=-1;
    bool found_slice = false;

    const rb::Cluster *theSlice = &(*fmSlice.at(i_track));

    for (size_t i_slice=0; i_slice<slices.size(); ++i_slice) {
      if (theSlice == slices.at(i_slice)){
      	slice = (float)i_slice; found_slice = true;
      	++tracks_per_slice.at(i_slice);
	//do the linear fit on the slice
      	if(containment != 4) slice_containment.at(i_slice) = 1;
      	break;
      }
    }

    if (!found_slice){
      slice = (float)slices.size();
      slices.push_back(theSlice);
      tracks_per_slice.push_back(1);
      if(containment!=4)
        slice_containment.push_back(1);
      else
        slice_containment.push_back(4);
    }

    // another loop to find slice id from vector slices
    slice = -1;
    for(size_t i_slice = 0; i_slice < slces.size(); i_slice++)
    {
    	if(slces.at(i_slice).get() == theSlice)
	{
	    slice = (float)i_slice;
	    break;
	}
    }


    if (!theTrack.Is3D()) continue;

    art::PtrVector<rb::CellHit> trackHits = theTrack.AllCells();
    std::set< std::pair<rb::CellHit,rb::RecoHit>,
      bool(*)(std::pair<rb::CellHit,rb::RecoHit>,
	      std::pair<rb::CellHit,rb::RecoHit>)>
      sortedTrackHits (pairHitComp);
    std::set< std::pair<rb::CellHit,rb::RecoHit>,
      bool(*)(std::pair<rb::CellHit,rb::RecoHit>,
	      std::pair<rb::CellHit,rb::RecoHit>)>
      sortedTrackHitsX (pairHitComp);
    std::set< std::pair<rb::CellHit,rb::RecoHit>,
      bool(*)(std::pair<rb::CellHit,rb::RecoHit>,
	      std::pair<rb::CellHit,rb::RecoHit>)>
      sortedTrackHitsY (pairHitComp);

    std::vector<double> x_hit;
    std::vector<double> y_hit;
    std::vector<double> zy_hit;
    std::vector<double> zx_hit;

    for (size_t i_hit=0; i_hit < trackHits.size(); ++i_hit) {
      art::Ptr<rb::CellHit> theHit = trackHits.at(i_hit);
      if (!theHit->GoodTiming()) continue;
      rb::RecoHit theRecoHit = theTrack.RecoHit(theHit);
      if (!theRecoHit.IsCalibrated()) continue;
      sortedTrackHits.insert(std::make_pair(*theHit,theRecoHit));

      unsigned short iC  = theHit->Cell();
      unsigned short iP  = theHit->Plane();
      geo::View_t   view = theHit->View();
      if(view == geo::kY){
        x_hit.push_back(iC);
        zx_hit.push_back(iP);
        sortedTrackHitsY.insert(std::make_pair(*theHit,theRecoHit));
      }
      else if(view == geo::kX){
        y_hit.push_back(iC);
        zy_hit.push_back(iP);
        sortedTrackHitsX.insert(std::make_pair(*theHit,theRecoHit));
      }
    }

    double fitpar[3];
    LinFit(x_hit, zx_hit, fitpar);
    double r2x = fitpar[2];
    LinFit(y_hit, zy_hit, fitpar);
    double r2y = fitpar[2];
    unsigned nxhit = x_hit.size();
    unsigned nyhit = y_hit.size();

    //////////////////////////////////////////////////////////////////////////
    //cp linear fit

   Double_t chi2X_fit = 1e6;
   Double_t chi2Y_fit = 1e6;
   const double wgt = 1;
   Double_t  trackX_fromfit = 0;
   Double_t trackY_fromfit = 0;

   Double_t  trackX_fromfit_slice = 0;
   Double_t trackY_fromfit_slice = 0;

   double  mvalX, cvalX,  mvalY, cvalY;
   TH1F *h1 = new TH1F("h1","xline",50,0,600);
   TH1F *h2 = new TH1F("h2","yline",50,0,600);

   //for (size_t i = 0; i<slices.size(); ++i) {
      std::vector<double> x_hit_fit;
      std::vector<double> y_hit_fit;
      std::vector<double> zy_hit_fit;
      std::vector<double> zx_hit_fit;
      
      std::vector<double> weight_x;
      std::vector<double> weight_y;


     //const rb::Cluster* slice = theSlice;

     art::PtrVector<rb::CellHit> xcells = theSlice->XCells();
     art::PtrVector<rb::CellHit> ycells = theSlice->YCells();

     for(size_t y_ind =0; y_ind < ycells.size(); y_ind++)
       {
	 art::Ptr<rb::CellHit> theHity = ycells.at(y_ind);
	
	 unsigned short iC  = theHity->Cell();
	 unsigned short iP  = theHity->Plane();
	  y_hit_fit.push_back(iC);
	  zy_hit_fit.push_back(iP);
	  weight_y.push_back(wgt);
	  trackY_fromfit += y_ind;
	  h1->Fill(y_hit_fit.at(y_ind));
	}

      for(size_t x_ind = 0; x_ind < xcells.size(); x_ind++)
	{
	  art::Ptr<rb::CellHit> theHitx = xcells.at(x_ind);

      	  unsigned short iC  = theHitx->Cell();
	  unsigned short iP  = theHitx->Plane();
	  x_hit_fit.push_back(iC);
	  zx_hit_fit.push_back(iP);
	  weight_x.push_back(wgt);
	  trackX_fromfit += x_ind;
	  h2->Fill(x_hit_fit.at(x_ind));
	}
      trackY_fromfit_slice += trackY_fromfit;
      trackX_fromfit_slice += trackX_fromfit;
    //}
       if(x_hit_fit.size()>=2){
      double chi2X_aux = util::LinFit(zx_hit_fit, x_hit_fit, weight_x, mvalX, cvalX);
      if(chi2X_aux==chi2X_aux)chi2X_fit = chi2X_aux;
      }
       if(y_hit_fit.size()>=2){
      double chi2Y_aux = util::LinFit(zy_hit_fit, y_hit_fit, weight_y, mvalY, cvalY);
      if(chi2Y_aux==chi2Y_aux)chi2Y_fit = chi2Y_aux;
      }

 ///////////////////////////////////////////////////////////
    std::vector<rb::RecoHit> outliers, outliersX, outliersY;
    double llr, P_up, P_dn, chi2, slope;
    double llrX, P_upX, P_dnX, chi2X, slopeX;
    double llrY, P_upY, P_dnY, chi2Y, slopeY;

    llr = getLLR(sortedTrackHits, outliers, P_up, P_dn, chi2, slope);
    llrX = getLLR(sortedTrackHitsX, outliersX, P_upX, P_dnX, chi2X, slopeX);
    llrY = getLLR(sortedTrackHitsY, outliersY, P_upY, P_dnY, chi2Y, slopeY);

    rb::RecoHit start;
    rb::RecoHit end;
    if (sortedTrackHits.size() >= 1) start = (sortedTrackHits.begin()->second);
    else { start = rb::RecoHit(); end = rb::RecoHit(); }
    if (sortedTrackHits.size() >= 2) end = (--sortedTrackHits.end())->second;
    else end = start;
    float dist = theTrack.TotalLength();

    // Average values
    float tDirX=0, tDirY=0, tDirZ=0, tEleAngle, trackTotalE=0;
    float avgTX=0, avgTY=0, tCosTheta, tAzimuth;

    for (size_t i_hit=0; i_hit < trackHits.size(); ++i_hit) {
      art::Ptr<rb::CellHit> theHit = trackHits.at(i_hit);
      rb::RecoHit theRecoHit = theTrack.RecoHit(theHit);

      TVector3 theDir = theTrack.InterpolateDir(theRecoHit.Z());

      tDirX += theDir.x()/trackHits.size();
      tDirY += theDir.y()/trackHits.size();
      tDirZ += theDir.z()/trackHits.size();

      if(theRecoHit.IsCalibrated()) trackTotalE += theRecoHit.GeV();

      if (theHit->View() == geo::kX)
        avgTX += theRecoHit.T()/sortedTrackHitsX.size();
      else if (theHit->View() == geo::kY)
        avgTY += theRecoHit.T()/sortedTrackHitsY.size();
    }

    tEleAngle = atan(tDirY/sqrt(tDirX*tDirX + tDirZ*tDirZ));
    tCosTheta = fUtilities.CosTheta(theTrack.Dir().Y(), theTrack.Dir().Mag());
    tAzimuth = fUtilities.Azimuth(theTrack.Dir().X(), theTrack.Dir().Z(), theTrack.Dir().Mag());

    TVector3 start_trk = theTrack.Start();
    TVector3 end_trk = theTrack.Stop();
    if(start_trk.Y() < end_trk.Y()) std::swap(start_trk, end_trk);
    const TVector3 dir_trk = (start_trk-end_trk).Unit();
    const double tzen_trk = acos(dir_trk.Y()) * 180 / M_PI;
    double tazi_trk = atan2(-dir_trk.X(), dir_trk.Z()) * 180 / M_PI + 332.066111;
    if(tazi_trk < 0) tazi_trk += 360;
    if(tazi_trk > 360) tazi_trk -= 360;

    const double tdotSun = AnglesToVector(zen, azi).Dot(AnglesToVector(tzen_trk, tazi_trk));

    float track_entries[55] =
      {
	(float)e.id().run(), (float)e.id().subRun(), (float)e.id().event(),
	slice, (float)i_track, (float)trackHits.size(), (float)sortedTrackHits.size(),
	(float)outliers.size(), (float)P_up, (float)P_dn, (float)llr,
        (float)chi2, (float)slope, (float)llrX, (float)chi2X, (float)slopeX,
        (float)llrY, (float)chi2Y, (float)slopeY, (float)r2x, (float)r2y,
        start.X(), start.Y(), start.Z(), start.T(), end.X(), end.Y(),
        end.Z(), end.T(), (float)nxhit,(float)nyhit, dist, tDirX, tDirY,
        tDirZ, tEleAngle, trackTotalE, containment, avgTX, avgTY,
	(float)raw_daqheader->TotalMicroSlices(), (float)zen, (float)azi,
	tCosTheta, tAzimuth, (float)tdotSun, (float)tzen_trk, (float)tazi_trk,
	(float)event_time, (float)chi2X_fit, (float)chi2Y_fit,(float)trackX_fromfit, (float)trackY_fromfit, (float)trackX_fromfit_slice, (float)trackY_fromfit_slice
      };
    ntp_track->Fill(track_entries);

    /*****************************************************************/
    /*  Matt Strait added this for the LIGO GW coincidence analysis  */
    /*****************************************************************/

    // Values from Cristiana on 2018-01-26
    const bool failscut =
         chi2 > 1.5
      || dist < 800.0 // Update from Cristiana on 2018-01-30
      || fabs(start.X() - end.X()) <  5*3.97 // 5 cell widths in cm
      || fabs(start.Y() - end.Y()) < 10*3.97
      || fabs(start.Z() - end.Z()) <  5*6.654 // 10 plane depths in cm
      || r2x < 0.99
      || r2y < 0.99
      || nxhit < 15
      || nyhit < 15
      || nxhit+nyhit < 60
      || llr < 3.0
      || slope < 0.0
      || slope > 2.0;

    if(!failscut){
      printf("Upmu: %d %d %d %.0f\t(%7.1f %7.1f %7.1f)\t(%7.1f %7.1f %7.1f)\n", e.run(), e.subRun(), e.event(), slice, start.X(), start.Y(), start.Z(), end.X(), end.Y(), end.Z());
      output->push_back(theTrack);
    }

    /*****************************************************************/
    /*    End section added for the LIGO GW coincidence analysis     */
    /*****************************************************************/
  } //end loop track
  e.put(std::move(output));
}

// Creates unit vector pointing to a point on unit sphere with theta=zen, phi = azi
TVector3 upmuana::MattUpMuAnalysis::AnglesToVector(double zen, double azi) const
{
  zen *= M_PI / 180;
  azi *= M_PI / 180;
  // Doesn't matter which axis is which
  return TVector3(sin(zen)*cos(azi), sin(zen)*sin(azi), cos(zen));
}

float upmuana::MattUpMuAnalysis::containmentType(rb::Track theTrack)
{
  TVector3 start = theTrack.Start();
  TVector3 stop  = theTrack.Stop();

  if (start.y() > stop.y()) // assume upward-going, so switch
    { TVector3 hold = start; start = stop; stop = hold; }

  if (start.y() < fMinY || start.x() < fMinX || start.x() > fMaxX ||
      start.z() < fMinZ || start.z() > fMaxZ) { // entered from outside
    if (stop.y() > fMaxY ||
        stop.x() < fMinX ||
        stop.x() > fMaxX ||
        stop.z() < fMinZ ||
        stop.z() > fMaxZ) { return 1; }//through-going
    return 2; // stopped
  } // started in detector
  if (stop.y() > fMaxY || stop.x() < fMinX || stop.x() > fMaxX ||
      stop.z() < fMinZ || stop.z() > fMaxZ)
    { return 3; }// in-produced

  return 4; // fully contained

}

double upmuana::MattUpMuAnalysis::getLLR(std::set< std::pair<rb::CellHit,rb::RecoHit>,
				    bool(*)(std::pair<rb::CellHit,rb::RecoHit>,
					    std::pair<rb::CellHit,rb::RecoHit>)>
				    sortedHits,
                                    std::vector<rb::RecoHit> &outliers,
                                    double &P_up, double &P_dn, double &chi2,
                                    double &slope)
{
  if (sortedHits.size() < 2) {
    P_up = 1e-30;
    P_dn = 1e-30;
    chi2 = 1e6;
    slope = -1e6;
    return 0;
  }

  std::vector<double> measuredTimes;
  std::vector<double> expectedTimes;
  std::vector<double> wgts;

  // first hit
  measuredTimes.push_back(0.0);
  expectedTimes.push_back(0.0);
  double err = getErr(sortedHits.begin()->first.PE());
  wgts.push_back(1.0/err/err);

  double  startX = sortedHits.begin()->second.X(),
    startY = sortedHits.begin()->second.Y(),
    startZ = sortedHits.begin()->second.Z(),
    startT = sortedHits.begin()->second.T();

  std::vector<rb::RecoHit> in;

  for (auto i_hit =  ++sortedHits.begin();
       i_hit != sortedHits.end();
       ++i_hit) {
    in.push_back(i_hit->second);
    measuredTimes.push_back(i_hit->second.T() - startT);
    double dist = getDist(i_hit->second.X(), i_hit->second.Y(),
                          i_hit->second.Z(), startX, startY, startZ);
    expectedTimes.push_back(dist/29.97);
    err = getErr(i_hit->first.PE());
    wgts.push_back(1.0/err/err);
  }

  LLR(expectedTimes, measuredTimes, wgts, slope, chi2, P_up, P_dn,
      in, outliers);

  return log(P_up/P_dn);
}

void upmuana::MattUpMuAnalysis::LLR(std::vector<double>& eT,
                               std::vector<double>& mT,
                               std::vector<double>& wgts, double& slope,
                               double& chi2, double& P_up, double& P_dn,
                               std::vector<rb::RecoHit>& in,
                               std::vector<rb::RecoHit>& outliers)
{
  // eT - x, mT - y
  size_t n = mT.size();

  // y = a + bx
  double a, b;
  if (eT.size() >= 2)
    chi2 = util::LinFit(eT, mT, wgts, b, a);
  else {
    chi2 = 1e6;
    slope = -1e6;
    P_up = 1e-30;
    P_dn = 1e-30;
    return;
  }

  size_t totAllowOutlier = n/10;
  size_t nOutliers = 0;
  std::vector<double> x_filt, y_filt, w_filt;
  for (size_t i=0; i < n; i++) {
    double y_fit = a + b*eT.at(i);
    if ( fabs(mT.at(i) - y_fit) < 5*(1.0/sqrt(wgts.at(i)))
         || nOutliers>totAllowOutlier)
      {
        x_filt.push_back(eT.at(i));
        y_filt.push_back(mT.at(i));
        w_filt.push_back(wgts.at(i));
      }
    else {
      ++nOutliers;
      outliers.push_back(in[i]);
    }
  }

  if (x_filt.size() >= 2)
    chi2 = util::LinFit(x_filt, y_filt, w_filt, b, a);
  else {
    chi2 = 1e6;
    slope = -1e6;
    P_up = 1e-30;
    P_dn = 1e-30;
  }
  n = x_filt.size();
  if (n < 5) {
    slope = 0;
    chi2  = 999;
    P_up  = 1e-30;
    P_dn  = 1;
    return;
  }

  slope = b;
  chi2 /= (double)(n-2);

  double  one_over_ee = 0.0,
    x_over_ee   = 0.0,
    y_over_ee   = 0.0;

  for (size_t i=0; i<n; ++i) {
    one_over_ee += w_filt.at(i);
    x_over_ee   += x_filt.at(i)*w_filt.at(i);
    y_over_ee   += y_filt.at(i)*w_filt.at(i);
  }

  double up_inter = (y_over_ee-x_over_ee)/one_over_ee;
  double dn_inter = (y_over_ee+x_over_ee)/one_over_ee;

  double up_chi2 = 0.0, dn_chi2 = 0.0;
  for (size_t i=0; i<n; ++i) {
    double e = 1.0/sqrt(w_filt.at(i));
    double up_expected = up_inter + x_filt.at(i);
    double dn_expected = dn_inter - x_filt.at(i);
    up_chi2 += pow((y_filt.at(i)-up_expected) / e, 2.0);
    dn_chi2 += pow((y_filt.at(i)-dn_expected) / e, 2.0);
  }

  double prob_up = boost::math::gamma_q((double)(n-2)/2.0,up_chi2/2.0),
    prob_dn = boost::math::gamma_q((double)(n-2)/2.0,dn_chi2/2.0);

  if (prob_up < 1e-30) prob_up = 1e-30;
  if (prob_dn < 1e-30) prob_dn = 1e-30;

  P_up = prob_up;
  P_dn = prob_dn;
}

void upmuana::MattUpMuAnalysis::LinFit(const std::vector<double>& x_val, const std::vector<double>& y_val, double *fitpar){

  const auto n    = x_val.size();
  const auto x  = (std::accumulate(x_val.begin(), x_val.end(), 0.0))/n;                       // <x>
  const auto y  = (std::accumulate(y_val.begin(), y_val.end(), 0.0))/n;                       // <y>
  const auto xx = (std::inner_product(x_val.begin(), x_val.end(), x_val.begin(), 0.0))/n;     // <xx>
  const auto yy = (std::inner_product(y_val.begin(), y_val.end(), y_val.begin(), 0.0))/n;     // <yy>
  const auto xy = (std::inner_product(x_val.begin(), x_val.end(), y_val.begin(), 0.0))/n;     // <xy>

  const auto b    = (xy - x*y)/(xx - x*x);                                                    // slope
  const auto a    = y - b*x;                                                                  // intercept
  const auto r    = (xy - x*y)/sqrt((xx - x*x)*(yy - y*y));                                   // Rxy - coeffcient of determination
  fitpar[0] = a;
  fitpar[1] = b;
  fitpar[2] = r*r;
}

void upmuana::MattUpMuAnalysis::LinFit(const std::vector<double>& x_val, const std::vector<double>& y_val, const std::vector<double>& y_err, double *fitpar){

  int n    = x_val.size();
  double x = 0;
  double y = 0;
  double xx = 0;
  double yy = 0;
  double xy = 0;
  double ee = 0;

  for ( int i=0; i<n; i++ ){
    x = x + x_val[i]/y_err[i]/y_err[i];
    y = y + y_val[i]/y_err[i]/y_err[i];

    xx = xx + x_val[i]*x_val[i]/y_err[i]/y_err[i];
    yy = yy + y_val[i]*y_val[i]/y_err[i]/y_err[i];
    xy = xy + x_val[i]*y_val[i]/y_err[i]/y_err[i];
    ee = ee + 1./y_err[i]/y_err[i];
  }

  const auto b    = (xy*ee - x*y)/(xx*ee - x*x);            // slope
  const auto a    = (xy - b*xx)/x;                          // intercept
  const auto r    = (xy - x*y)/sqrt((xx - x*x)*(yy - y*y)); // Rxy - coeffcient of determination
  fitpar[0] = a;
  fitpar[1] = b;
  fitpar[2] = r*r;
}

DEFINE_ART_MODULE(upmuana::MattUpMuAnalysis)
