////////////////////////////////////////////////////////////////////////
// Class:       RecoDump
// Plugin Type: analyzer (Unknown Unknown)
// File:        RecoDump.cc
//
// Tue 12 Jul 22 Alex Wilkinson
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

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
// #include "nusimdata/SimulationBase/MCTruth.h"
// #include "nusimdata/SimulationBase/MCParticle.h"
// #include "nusimdata/SimulationBase/MCNeutrino.h"
// #include "lardataobj/AnalysisBase/ParticleID.h"
// #include "lardataobj/RecoBase/Track.h"
// #include "lardataobj/RecoBase/PFParticle.h"
// #include "lardataobj/RecoBase/Shower.h"
// #include "lardataobj/AnalysisBase/Calorimetry.h"
// #include "lardataobj/RecoBase/Wire.h"
#include "dunereco/CVN/func/Result.h"
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include <string>

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace extrapolation {
  class RecoDump;
}


class extrapolation::RecoDump : public art::EDAnalyzer {
public:
  explicit RecoDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RecoDump(RecoDump const&) = delete;
  RecoDump(RecoDump&&) = delete;
  RecoDump& operator=(RecoDump const&) = delete;
  RecoDump& operator=(RecoDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My function.
  void reset();

private:
  const geo::GeometryCore* fGeom;

  // Labels from fcl
  std::string fTrueCVNResultsLabel;
  std::string fNetworkCVNResultsLabel;
  std::string fTrueNumuEResultsLabel;
  std::string fNetworkNumuEResultsLabel;
  std::string fTrueNCEResultsLabel;
  std::string fNetworkNCEResultsLabel;
  std::string fEventIDSEDLabel;

  TTree*  fTreeReco;
  int     fEventID;
  int     fRun;
  int     fSubRun;
  int     fEventNum;
  // Flavour scores
  float fTrueNumuScore;
  float fNetworkNumuScore;
  float fTrueNueScore;
  float fNetworkNueScore;
  float fTrueNCScore;
  float fNetworkNCScore;
  float fTrueNutauScore;
  float fNetworkNutauScore;
  // Antinuetrino score
  float fTrueAntiNuScore;
  float fNetworkAntiNuScore;
  // Nu energy reco information
  double fTrueNumuNuE;
  double fNetworkNumuNuE;
  double fTrueNumuHadE;
  double fNetworkNumuHadE;
  double fTrueNumuLepE;
  double fNetworkNumuLepE;
  double fTrueNCNuE;
  double fNetworkNCNuE;
  double fTrueNCHadE;
  double fNetworkNCHadE;
  double fTrueNCLepE;
  double fNetworkNCLepE;
};


extrapolation::RecoDump::RecoDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTrueCVNResultsLabel      (p.get<std::string>("TrueCVNResultsLabel")),
    fNetworkCVNResultsLabel   (p.get<std::string>("NetworkCVNResultsLabel")),
    fTrueNumuEResultsLabel    (p.get<std::string>("TrueNumuEResultsLabel")),
    fNetworkNumuEResultsLabel (p.get<std::string>("NetworkNumuEResultsLabel")),
    fTrueNCEResultsLabel      (p.get<std::string>("TrueNCEResultsLabel")),
    fNetworkNCEResultsLabel   (p.get<std::string>("NetworkNCEResultsLabel")),
    fEventIDSEDLabel          (p.get<std::string>("EventIDSEDLabel"))
{
  consumes<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);

  consumes<std::vector<cvn::Result>>(fTrueCVNResultsLabel);
  consumes<std::vector<cvn::Result>>(fNetworkCVNResultsLabel);

  consumes<dune::EnergyRecoOutput>(fTrueNumuEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNetworkNumuEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fTrueNCEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNetworkNCEResultsLabel);

  art::ServiceHandle<art::TFileService> tfs;

  fTreeReco = tfs->make<TTree>("FDReco", "FDReco");
  fTreeReco->Branch("eventID", &fEventID, "eventid/I");
  fTreeReco->Branch("Run", &fRun, "run/I");
  fTreeReco->Branch("SubRun", &fSubRun, "subrun/I");
  fTreeReco->Branch("EventNum", &fEventNum, "eventnum/I");
  // Flavour scores
  fTreeReco->Branch("TrueNumuScore", &fTrueNumuScore, "truenumuscore/F");
  fTreeReco->Branch("NetworkNumuScore", &fNetworkNumuScore, "networknumuscore/F");
  fTreeReco->Branch("TrueNueScore", &fTrueNueScore, "truenuescore/F");
  fTreeReco->Branch("NetworkNueScore", &fNetworkNueScore, "networknuescore/F");
  fTreeReco->Branch("TrueNCScore", &fTrueNCScore, "truencscore/F");
  fTreeReco->Branch("NetworkNCScore", &fNetworkNCScore, "networkncscore/F");
  fTreeReco->Branch("TrueNutauScore", &fTrueNutauScore, "truenutauscore/F");
  fTreeReco->Branch("NetworkNutauScore", &fNetworkNutauScore, "networknutauscore/F");
  // Antineutrino score
  fTreeReco->Branch("TrueAntiNuScore", &fTrueAntiNuScore, "trueantinuscore/F");
  fTreeReco->Branch("NetworkAntiNuScore", &fNetworkAntiNuScore, "networkantinuscore/F");
  // Nu energy reco information
  fTreeReco->Branch("TrueNumuNuE", &fTrueNumuNuE, "truenumunue/D");
  fTreeReco->Branch("NetworkNumuNuE", &fNetworkNumuNuE, "networknumunue/D");
  fTreeReco->Branch("TrueNumuHadE", &fTrueNumuHadE, "truenumuhade/D");
  fTreeReco->Branch("NetworkNumuHadE", &fNetworkNumuHadE, "networknumuhade/D");
  fTreeReco->Branch("TrueNumuLepE", &fTrueNumuLepE, "truenumulepe/D");
  fTreeReco->Branch("NetworkNumuLepE", &fNetworkNumuLepE, "networknumulepe/D");
  fTreeReco->Branch("TrueNCNuE", &fTrueNCNuE, "truencnue/D");
  fTreeReco->Branch("NetworkNCNuE", &fNetworkNCNuE, "networkncnue/D");
  fTreeReco->Branch("TrueNCHadE", &fTrueNCHadE, "truenchade/D");
  fTreeReco->Branch("NetworkNCHadE", &fNetworkNCHadE, "networknchade/D");
  fTreeReco->Branch("TrueNCLepE", &fTrueNCLepE, "truenclepe/D");
  fTreeReco->Branch("NetworkNCLepE", &fNetworkNCLepE, "networknclepe/D");
}

void extrapolation::RecoDump::analyze(art::Event const& e)
{
  this->reset();

  // Get SED containing eventID
  const auto eventIDSEDs = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);
  fEventID = eventIDSEDs->at(0).TrackID();

  // Get event id info
  fRun = e.id().run();
  fSubRun = e.id().subRun();
  fEventNum = e.id().event();

  // Get results from CVN
  const auto trueCVNResults = e.getValidHandle<std::vector<cvn::Result>>(fTrueCVNResultsLabel);
  const auto networkCVNResults = e.getValidHandle<std::vector<cvn::Result>>(fNetworkCVNResultsLabel);

  // Get flavour scores
  fTrueNumuScore = trueCVNResults->at(0).GetNumuProbability();
  fNetworkNumuScore = networkCVNResults->at(0).GetNumuProbability();
  fTrueNueScore = trueCVNResults->at(0).GetNueProbability();
  fNetworkNueScore = networkCVNResults->at(0).GetNueProbability();
  fTrueNCScore = trueCVNResults->at(0).GetNCProbability();
  fNetworkNCScore = networkCVNResults->at(0).GetNCProbability();
  fTrueNutauScore = trueCVNResults->at(0).GetNutauProbability();
  fNetworkNutauScore = networkCVNResults->at(0).GetNutauProbability();

  // Get antineutrino score
  fTrueAntiNuScore = trueCVNResults->at(0).GetIsAntineutrinoProbability();
  fNetworkAntiNuScore = networkCVNResults->at(0).GetIsAntineutrinoProbability();

  // Get nu reco information
  const auto trueNumuEOut = e.getValidHandle<dune::EnergyRecoOutput> (fTrueNumuEResultsLabel);
  const auto networkNumuEOut = e.getValidHandle<dune::EnergyRecoOutput> (fNetworkNumuEResultsLabel);
  const auto trueNCEOut = e.getValidHandle<dune::EnergyRecoOutput> (fTrueNCEResultsLabel);
  const auto networkNCEOut = e.getValidHandle<dune::EnergyRecoOutput> (fNetworkNCEResultsLabel);

  fTrueNumuNuE = trueNumuEOut->fNuLorentzVector.E();
  fNetworkNumuNuE = networkNumuEOut->fNuLorentzVector.E();
  fTrueNumuHadE = trueNumuEOut->fHadLorentzVector.E();
  fNetworkNumuHadE = networkNumuEOut->fHadLorentzVector.E();
  fTrueNumuLepE = trueNumuEOut->fLepLorentzVector.E();
  fNetworkNumuLepE = networkNumuEOut->fLepLorentzVector.E();

  fTrueNCNuE = trueNCEOut->fNuLorentzVector.E();
  fNetworkNCNuE = networkNCEOut->fNuLorentzVector.E();
  fTrueNCHadE = trueNCEOut->fHadLorentzVector.E();
  fNetworkNCHadE = networkNCEOut->fHadLorentzVector.E();
  fTrueNCLepE = trueNCEOut->fLepLorentzVector.E();
  fNetworkNCLepE = networkNCEOut->fLepLorentzVector.E();

  fTreeReco->Fill();
}

void extrapolation::RecoDump::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
}

void extrapolation::RecoDump::endJob()
{
}

void extrapolation::RecoDump::reset()
{
  fEventID = -99;
  fRun = -999;
  fSubRun = -999;
  fEventNum = -999;

  fTrueNumuScore = -999.0;
  fNetworkNumuScore = -999.0;
  fTrueNueScore = -999.0;
  fNetworkNueScore = -999.0;
  fTrueNCScore = -999.0;
  fNetworkNCScore = -999.0;
  fTrueNutauScore = -999.0;
  fNetworkNutauScore = -999.0;


  fTrueAntiNuScore = -999.0;
  fNetworkAntiNuScore = -999.0;

  fTrueNumuNuE = -1.0;
  fNetworkNumuNuE = -1.0;
  fTrueNumuHadE = -1.0;
  fNetworkNumuHadE = -1.0;
  fTrueNumuLepE = -1.0;
  fNetworkNumuLepE = -1.0;
  fTrueNCNuE = -1.0;
  fNetworkNCNuE = -1.0;
  fTrueNCHadE = -1.0;
  fNetworkNCHadE = -1.0;
  fTrueNCLepE = -1.0;
  fNetworkNCLepE = -1.0;
}

DEFINE_ART_MODULE(extrapolation::RecoDump)
