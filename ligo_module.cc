////////////////////////////////////////////////////////////////////////
/// \brief   This module is named ligo and looks at LIGO coincidences
/// \author  M. Strait
/// \date
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "RecoBase/Track.h"
#include "Calibrator/Calibrator.h"

#include "Simulation/ParticleNavigator.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCNeutrino.h"

#include <string>
#include <algorithm>

#include <signal.h>

#include "art/Framework/Core/EDAnalyzer.h"


/// Calibrating RawData to Produce CellHits
namespace ligo {

  class ligo : public art::EDAnalyzer {

  public:

    explicit ligo(fhicl::ParameterSet const& pset);
    virtual ~ligo();

    void analyze(const art::Event& evt);

  }; // class ligo
}


namespace ligo{

ligo::ligo(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
{ }

ligo::~ligo() { }

void ligo::analyze(const art::Event& evt)
{
  // I'm so sorry that I have to do this.  And, my goodness, doing
  // it in the constructor isn't sufficient.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<sim::Particle> > particles;
  evt.getByLabel("geantgen", particles);

  {
    static int NOvA = printf(
      "ntuple: totallength\n");
    NOvA = NOvA;
  }

  printf("ntuple: %.1f", 
                    (*particles)[0].Trajectory().TotalLength()
        );

  printf("\n");
}

DEFINE_ART_MODULE(ligo);

} // end namespace ligo
//////////////////////////////////////////////////////////////////////////
