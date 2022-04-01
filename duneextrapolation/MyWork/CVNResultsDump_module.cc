////////////////////////////////////////////////////////////////////////
// Class:       CVNResultsDump
// Plugin Type: analyzer (Unknown Unknown)
// File:        CVNResultsDump_module.cc
//
// Thur Mar 31 2022 Alex Wilkinson
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "dunereco/CVN/func/Result.h"

#include <string>

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace extrapolation {
  class CVNResultsDump;
}


class extrapolation::CVNResultsDump : public art::EDAnalyzer {
public:
  explicit CVNResultsDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CVNResultsDump(CVNResultsDump const&) = delete;
  CVNResultsDump(CVNResultsDump&&) = delete;
  CVNResultsDump& operator=(CVNResultsDump const&) = delete;
  CVNResultsDump& operator=(CVNResultsDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My function.
  void reset();

private:
  const geo::GeometryCore* fGeom;

  std::string fTrueCVNResultsLabel;
  std::string fNetworkCVNResultsLabel;
  std::string fNDCVNResultsLabel;
  
  TTree* fTreeCVNResults;
  int    fRun;
  int    fSubRun;
  int    fEventNum;
  // Flavour scores
  float  fTrueNumuScore;
  float  fNetworkNumuScore;
  float  fNDNumuScore;
  float  fTrueNueScore;
  float  fNetworkNueScore;
  float  fNDNueScore;
  float  fTrueNCScore;
  float  fNetworkNCScore;
  float  fNDNCScore;
  float  fTrueNutauScore;
  float  fNetworkNutauScore;
  float  fNDNutauScore;
  // Pred number of daughter particle types
  int    fTrueNumPions;
  int    fNetworkNumPions;
  int    fNDNumPions;
  int    fTrueNumProtons;
  int    fNetworkNumProtons;
  int    fNDNumProtons;
  int    fTrueNumPizeros;
  int    fNetworkNumPizeros;
  int    fNDNumPizeros;
  int    fTrueNumNeutrons;
  int    fNetworkNumNeutrons;
  int    fNDNumNeutrons;
  // Antinuetrino score
  float  fTrueAntiNuScore;
  float  fNetworkAntiNuScore;
  float  fNDAntiNuScore;
  // Interaction type score
  float  fTrueQEScore;
  float  fNetworkQEScore;
  float  fNDQEScore;
  float  fTrueResScore;
  float  fNetworkResScore;
  float  fNDResScore;
  float  fTrueDISScore;
  float  fNetworkDISScore;
  float  fNDDISScore;
  float  fTrueOtherScore;
  float  fNetworkOtherScore;
  float  fNDOtherScore;
};


extrapolation::CVNResultsDump::CVNResultsDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTrueCVNResultsLabel       (p.get<std::string> ("TrueCVNResultsLabel")),
    fNetworkCVNResultsLabel (p.get<std::string> ("NetworkCVNResultsLabel")),
    fNDCVNResultsLabel         (p.get<std::string> ("NDCVNResultsLabel"))
{
  consumes<std::vector<cvn::Result>>(fTrueCVNResultsLabel);
  consumes<std::vector<cvn::Result>>(fNetworkCVNResultsLabel);
  consumes<std::vector<cvn::Result>>(fNDCVNResultsLabel);

  art::ServiceHandle<art::TFileService> tfs;

  fTreeCVNResults = tfs->make<TTree>("CVNResults", "CVNResults");
  fTreeCVNResults->Branch("Run", &fRun, "run/I");
  fTreeCVNResults->Branch("SubRun", &fSubRun, "subrun/I");
  fTreeCVNResults->Branch("EventNum", &fEventNum, "eventnum/I");
  // Flavour scores
  fTreeCVNResults->Branch("TrueNumuScore", &fTrueNumuScore, "truenumuscore/F");
  fTreeCVNResults->Branch("NetworkNumuScore", &fNetworkNumuScore, "networknumuscore/F");
  fTreeCVNResults->Branch("NDNumuScore", &fNDNumuScore, "ndnumuscore/F");
  fTreeCVNResults->Branch("TrueNueScore", &fTrueNueScore, "truenuescore/F");
  fTreeCVNResults->Branch("NetworkNueScore", &fNetworkNueScore, "networknuescore/F");
  fTreeCVNResults->Branch("NDNueScore", &fNDNueScore, "ndnuescore/F");  
  fTreeCVNResults->Branch("TrueNCScore", &fTrueNCScore, "truencscore/F");
  fTreeCVNResults->Branch("NetworkNCScore", &fNetworkNCScore, "networkncscore/F");
  fTreeCVNResults->Branch("NDNCScore", &fNDNCScore, "ndncscore/F");
  fTreeCVNResults->Branch("TrueNutauScore", &fTrueNutauScore, "truenutauscore/F");
  fTreeCVNResults->Branch("NetworkNutauScore", &fNetworkNutauScore, "networknutauscore/F");
  fTreeCVNResults->Branch("NDNutauScore", &fNDNutauScore, "ndnutauscore/F");
  // Pred number of daughter particles
  fTreeCVNResults->Branch("TrueNumPions", &fTrueNumPions, "truenumpions/I");
  fTreeCVNResults->Branch("NetworkNumPions", &fNetworkNumPions, "networknumpions/I");
  fTreeCVNResults->Branch("NDNumPions", &fNDNumPions, "ndnumpions/I");
  fTreeCVNResults->Branch("TrueNumProtons", &fTrueNumProtons, "truenumprotons/I");
  fTreeCVNResults->Branch("NetworkNumProtons", &fNetworkNumProtons, "networknumprotons/I");
  fTreeCVNResults->Branch("NDNumProtons", &fNDNumProtons, "ndnumprotons/I");
  fTreeCVNResults->Branch("TrueNumPizeros", &fTrueNumPizeros, "truenumpizeros/I");
  fTreeCVNResults->Branch("NetworkNumPizeros", &fNetworkNumPizeros, "networknumpizeros/I");
  fTreeCVNResults->Branch("NDNumPizeros", &fNDNumPizeros, "ndnumpizeros/I");
  fTreeCVNResults->Branch("TrueNumNeutrons", &fTrueNumNeutrons, "truenumneutrons/I");
  fTreeCVNResults->Branch("NetworkNumNeutrons", &fNetworkNumNeutrons, "networknumneutrons/I");
  fTreeCVNResults->Branch("NDNumNeutrons", &fNDNumNeutrons, "ndnumneutrons/I");
  // Antineutrino score
  fTreeCVNResults->Branch("TrueAntiNuScore", &fTrueAntiNuScore, "trueantinuscore/F");
  fTreeCVNResults->Branch("NetworkAntiNuScore", &fNetworkAntiNuScore, "networkantinuscore/F");
  fTreeCVNResults->Branch("NDAntiNuScore", &fNDAntiNuScore, "ndantinuscore/F");
  // Interaction Type scores
  fTreeCVNResults->Branch("TrueQEScore", &fTrueQEScore, "trueqescore/F");
  fTreeCVNResults->Branch("NetworkQEScore", &fNetworkQEScore, "networkqescore/F");
  fTreeCVNResults->Branch("NDQEScore", &fNDQEScore, "ndqescore/F");
  fTreeCVNResults->Branch("TrueResScore", &fTrueResScore, "trueresscore/F");
  fTreeCVNResults->Branch("NetworkResScore", &fNetworkResScore, "networkresscore/F");
  fTreeCVNResults->Branch("NDResScore", &fNDResScore, "ndresscore/F");
  fTreeCVNResults->Branch("TrueDISScore", &fTrueDISScore, "truedisscore/F");
  fTreeCVNResults->Branch("NetworkDISScore", &fNetworkDISScore, "networkdisscore/F");
  fTreeCVNResults->Branch("NDDISScore", &fNDDISScore, "nddisscore/F");
  fTreeCVNResults->Branch("TrueOtherScore", &fTrueOtherScore, "trueotherscore/F");
  fTreeCVNResults->Branch("NetworkOtherScore", &fNetworkOtherScore, "networkotherscore/F");
  fTreeCVNResults->Branch("NDOtherScore", &fNDQEScore, "ndotherscore/F"); 
}

void extrapolation::CVNResultsDump::analyze(art::Event const& e)
{
  this->reset();

  // Get event id info
  fRun = e.id().run();
  fSubRun = e.id().subRun(); 
  fEventNum = e.id().event();

  const auto trueCVNResults = e.getValidHandle<std::vector<cvn::Result>> (fTrueCVNResultsLabel);
  const auto networkCVNResults = e.getValidHandle<std::vector<cvn::Result>> (fNetworkCVNResultsLabel);
  const auto NDCVNResults = e.getValidHandle<std::vector<cvn::Result>> (fNDCVNResultsLabel);
  
  // Get flavour scores
  fTrueNumuScore = trueCVNResults->at(0).GetNumuProbability();
  fNetworkNumuScore = networkCVNResults->at(0).GetNumuProbability();
  fNDNumuScore = NDCVNResults->at(0).GetNumuProbability();
  fTrueNueScore = trueCVNResults->at(0).GetNueProbability();
  fNetworkNueScore = networkCVNResults->at(0).GetNueProbability();
  fNDNueScore = NDCVNResults->at(0).GetNueProbability();
  fTrueNCScore = trueCVNResults->at(0).GetNCProbability();
  fNetworkNCScore = networkCVNResults->at(0).GetNCProbability();
  fNDNCScore = NDCVNResults->at(0).GetNCProbability();
  fTrueNutauScore = trueCVNResults->at(0).GetNutauProbability();
  fNetworkNutauScore = networkCVNResults->at(0).GetNutauProbability();
  fNDNutauScore = NDCVNResults->at(0).GetNutauProbability();
  
  // Get pred number of daughter particle types (3 means >2)
  fTrueNumPions = (int)trueCVNResults->at(0).PredictedPions();
  fNetworkNumPions = (int)networkCVNResults->at(0).PredictedPions();
  fNDNumPions = (int)NDCVNResults->at(0).PredictedPions();
  fTrueNumProtons = (int)trueCVNResults->at(0).PredictedProtons();
  fNetworkNumProtons = (int)networkCVNResults->at(0).PredictedProtons();
  fNDNumProtons = (int)NDCVNResults->at(0).PredictedProtons();
  fTrueNumPizeros = (int)trueCVNResults->at(0).PredictedPizeros();
  fNetworkNumPizeros = (int)networkCVNResults->at(0).PredictedPizeros();
  fNDNumPizeros = (int)NDCVNResults->at(0).PredictedPizeros();
  fTrueNumNeutrons = (int)trueCVNResults->at(0).PredictedNeutrons();
  fNetworkNumNeutrons = (int)networkCVNResults->at(0).PredictedNeutrons();
  fNDNumNeutrons = (int)NDCVNResults->at(0).PredictedNeutrons();

  // Get antineutrino score
  fTrueAntiNuScore = trueCVNResults->at(0).GetIsAntineutrinoProbability();
  fNetworkAntiNuScore = networkCVNResults->at(0).GetIsAntineutrinoProbability();
  fNDAntiNuScore = NDCVNResults->at(0).GetIsAntineutrinoProbability();

  // Get interaction type scores
  fTrueQEScore = trueCVNResults->at(0).GetQEProbability();
  fNetworkQEScore = networkCVNResults->at(0).GetQEProbability();
  fNDQEScore = NDCVNResults->at(0).GetQEProbability();
  fTrueResScore = trueCVNResults->at(0).GetResProbability();
  fNetworkResScore = networkCVNResults->at(0).GetResProbability();
  fNDResScore = NDCVNResults->at(0).GetResProbability();
  fTrueDISScore = trueCVNResults->at(0).GetDISProbability();
  fNetworkDISScore = networkCVNResults->at(0).GetDISProbability();
  fNDDISScore = NDCVNResults->at(0).GetDISProbability();
  fTrueOtherScore = trueCVNResults->at(0).GetOtherProbability();
  fNetworkOtherScore = networkCVNResults->at(0).GetOtherProbability();
  fNDOtherScore = NDCVNResults->at(0).GetOtherProbability();
  
  fTreeCVNResults->Fill();
}

void extrapolation::CVNResultsDump::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
}

void extrapolation::CVNResultsDump::endJob()
{
}

void extrapolation::CVNResultsDump::reset() 
{
  fRun = -999;
  fSubRun = -999;
  fEventNum = -999;

  fTrueNumuScore = -999.0;
  fNetworkNumuScore = -999.0;
  fNDNumuScore = -999.0;
  fTrueNueScore = -999.0;
  fNetworkNueScore = -999.0;
  fNDNueScore = -999.0;
  fTrueNCScore = -999.0;
  fNetworkNCScore = -999.0;
  fNDNCScore = -999.0;
  fTrueNutauScore = -999.0;
  fNetworkNutauScore = -999.0;
  fNDNutauScore = -999.0;

  fTrueNumPions = -999;
  fNetworkNumPions = -999;
  fNDNumPions = -999;
  fTrueNumProtons = -999;
  fNetworkNumProtons = -999;
  fNDNumProtons = -999;
  fTrueNumPizeros = -999;
  fNetworkNumPizeros = -999;
  fNDNumPizeros = -999;
  fTrueNumNeutrons = -999;
  fNetworkNumNeutrons = -999;
  fNDNumNeutrons = -999;

  fTrueAntiNuScore = -999.0;
  fNetworkAntiNuScore = -999.0;
  fNDAntiNuScore = -999.0;

  fTrueQEScore = -999.0;
  fNetworkQEScore = -999.0;
  fNDQEScore = -999.0;
  fTrueResScore = -999.0;
  fNetworkResScore = -999.0;
  fNDResScore = -999.0;
  fTrueDISScore = -999.0;
  fNetworkDISScore = -999.0;
  fNDDISScore = -999.0;
  fTrueOtherScore = -999.0;
  fNetworkOtherScore = -999.0;
  fNDOtherScore = -999.0;
}

DEFINE_ART_MODULE(extrapolation::CVNResultsDump)
