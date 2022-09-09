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
#include "dunereco/CVN/func/Result.h"
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/Hit.h"

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
  std::string fTrueNueEResultsLabel;
  std::string fNetworkNueEResultsLabel;
  std::string fTrueNCEResultsLabel;
  std::string fNetworkNCEResultsLabel;
  std::string fTrueHitsLabel;
  std::string fNetworkHitsLabel;
  std::string fEventIDSEDLabel;

  TTree* fTreeReco;
  int    fEventID;
  int    fRun;
  int    fSubRun;
  int    fEventNum;
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
  float fTrueNumuNuE;
  float fNetworkNumuNuE;
  float fTrueNumuHadE;
  float fNetworkNumuHadE;
  float fTrueNumuLepE;
  float fNetworkNumuLepE;
  float fTrueNueNuE;
  float fNetworkNueNuE;
  float fTrueNueHadE;
  float fNetworkNueHadE;
  float fTrueNueLepE;
  float fNetworkNueLepE;
  float fTrueNCNuE;
  float fNetworkNCNuE;
  float fTrueNCHadE;
  float fNetworkNCHadE;
  float fTrueNCLepE;
  float fNetworkNCLepE;
  // Hit information
  int fTrueNumHits;
  int fNetworkNumHits;

  int fNumOORErrs;
};


extrapolation::RecoDump::RecoDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTrueCVNResultsLabel      (p.get<std::string>("TrueCVNResultsLabel")),
    fNetworkCVNResultsLabel   (p.get<std::string>("NetworkCVNResultsLabel")),
    fTrueNumuEResultsLabel    (p.get<std::string>("TrueNumuEResultsLabel")),
    fNetworkNumuEResultsLabel (p.get<std::string>("NetworkNumuEResultsLabel")),
    fTrueNueEResultsLabel     (p.get<std::string>("TrueNueEResultsLabel")),
    fNetworkNueEResultsLabel  (p.get<std::string>("NetworkNueEResultsLabel")),
    fTrueNCEResultsLabel      (p.get<std::string>("TrueNCEResultsLabel")),
    fNetworkNCEResultsLabel   (p.get<std::string>("NetworkNCEResultsLabel")),
    fTrueHitsLabel            (p.get<std::string>("TrueHitsLabel")),
    fNetworkHitsLabel         (p.get<std::string>("NetworkHitsLabel")),
    fEventIDSEDLabel          (p.get<std::string>("EventIDSEDLabel"))
{
  consumes<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);

  consumes<std::vector<cvn::Result>>(fTrueCVNResultsLabel);
  consumes<std::vector<cvn::Result>>(fNetworkCVNResultsLabel);

  consumes<dune::EnergyRecoOutput>(fTrueNumuEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNetworkNumuEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fTrueNueEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNetworkNueEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fTrueNCEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNetworkNCEResultsLabel);

  consumes<std::vector<recob::Hit>>(fTrueHitsLabel);
  consumes<std::vector<recob::Hit>>(fNetworkHitsLabel);

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
  fTreeReco->Branch("TrueNumuNuE", &fTrueNumuNuE, "truenumunue/F");
  fTreeReco->Branch("NetworkNumuNuE", &fNetworkNumuNuE, "networknumunue/F");
  fTreeReco->Branch("TrueNumuHadE", &fTrueNumuHadE, "truenumuhade/F");
  fTreeReco->Branch("NetworkNumuHadE", &fNetworkNumuHadE, "networknumuhade/F");
  fTreeReco->Branch("TrueNumuLepE", &fTrueNumuLepE, "truenumulepe/F");
  fTreeReco->Branch("NetworkNumuLepE", &fNetworkNumuLepE, "networknumulepe/F");
  fTreeReco->Branch("TrueNueNuE", &fTrueNueNuE, "truenuenue/F");
  fTreeReco->Branch("NetworkNueNuE", &fNetworkNueNuE, "networknuenue/F");
  fTreeReco->Branch("TrueNueHadE", &fTrueNueHadE, "truenuehade/F");
  fTreeReco->Branch("NetworkNueHadE", &fNetworkNueHadE, "networknuehade/F");
  fTreeReco->Branch("TrueNueLepE", &fTrueNueLepE, "truenuelepe/F");
  fTreeReco->Branch("NetworkNueLepE", &fNetworkNueLepE, "networknuelepe/F");
  fTreeReco->Branch("TrueNCNuE", &fTrueNCNuE, "truencnue/F");
  fTreeReco->Branch("NetworkNCNuE", &fNetworkNCNuE, "networkncnue/F");
  fTreeReco->Branch("TrueNCHadE", &fTrueNCHadE, "truenchade/F");
  fTreeReco->Branch("NetworkNCHadE", &fNetworkNCHadE, "networknchade/F");
  fTreeReco->Branch("TrueNCLepE", &fTrueNCLepE, "truenclepe/F");
  fTreeReco->Branch("NetworkNCLepE", &fNetworkNCLepE, "networknclepe/F");
  // Hit information
  fTreeReco->Branch("TrueNumHits", &fTrueNumHits, "truenumhits/I");
  fTreeReco->Branch("NetworkNumHits", &fNetworkNumHits, "networknumhits/I");

}

void extrapolation::RecoDump::analyze(art::Event const& e)
{
  this->reset();

  try {
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
    const auto trueNumuEOut = e.getValidHandle<dune::EnergyRecoOutput>(fTrueNumuEResultsLabel);
    const auto networkNumuEOut = e.getValidHandle<dune::EnergyRecoOutput>(fNetworkNumuEResultsLabel);
    const auto trueNueEOut = e.getValidHandle<dune::EnergyRecoOutput>(fTrueNueEResultsLabel);
    const auto networkNueEOut = e.getValidHandle<dune::EnergyRecoOutput>(fNetworkNueEResultsLabel);
    const auto trueNCEOut = e.getValidHandle<dune::EnergyRecoOutput>(fTrueNCEResultsLabel);
    const auto networkNCEOut = e.getValidHandle<dune::EnergyRecoOutput>(fNetworkNCEResultsLabel);

    fTrueNumuNuE = (float)trueNumuEOut->fNuLorentzVector.E();
    fNetworkNumuNuE = (float)networkNumuEOut->fNuLorentzVector.E();
    fTrueNumuHadE = (float)trueNumuEOut->fHadLorentzVector.E();
    fNetworkNumuHadE = (float)networkNumuEOut->fHadLorentzVector.E();
    fTrueNumuLepE = (float)trueNumuEOut->fLepLorentzVector.E();
    fNetworkNumuLepE = (float)networkNumuEOut->fLepLorentzVector.E();

    fTrueNueNuE = (float)trueNueEOut->fNuLorentzVector.E();
    fNetworkNueNuE = (float)networkNueEOut->fNuLorentzVector.E();
    fTrueNueHadE = (float)trueNueEOut->fHadLorentzVector.E();
    fNetworkNueHadE = (float)networkNueEOut->fHadLorentzVector.E();
    fTrueNueLepE = (float)trueNueEOut->fLepLorentzVector.E();
    fNetworkNueLepE = (float)networkNueEOut->fLepLorentzVector.E();

    fTrueNCNuE = (float)trueNCEOut->fNuLorentzVector.E();
    fNetworkNCNuE = (float)networkNCEOut->fNuLorentzVector.E();
    fTrueNCHadE = (float)trueNCEOut->fHadLorentzVector.E();
    fNetworkNCHadE = (float)networkNCEOut->fHadLorentzVector.E();
    fTrueNCLepE = (float)trueNCEOut->fLepLorentzVector.E();
    fNetworkNCLepE = (float)networkNCEOut->fLepLorentzVector.E();

    // Get hit information
    const auto trueHits = e.getValidHandle<std::vector<recob::Hit>>(fTrueHitsLabel);
    const auto networkHits = e.getValidHandle<std::vector<recob::Hit>>(fNetworkHitsLabel);

    fTrueNumHits = (int)trueHits->size();
    fNetworkNumHits = (int)networkHits->size();
  }
  catch (const std::out_of_range& err) {
    fNumOORErrs++;

    this->reset();
    fTreeReco->Fill();
  }
}

void extrapolation::RecoDump::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  fNumOORErrs = 0;
}

void extrapolation::RecoDump::endJob()
{
  std::cout << "There were " << fNumOORErrs << 
               " out of range errors when accessing data products\n";
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
  fTrueNueNuE = -1.0;
  fNetworkNueNuE = -1.0;
  fTrueNueHadE = -1.0;
  fNetworkNueHadE = -1.0;
  fTrueNueLepE = -1.0;
  fNetworkNueLepE = -1.0;
  fTrueNCNuE = -1.0;
  fNetworkNCNuE = -1.0;
  fTrueNCHadE = -1.0;
  fNetworkNCHadE = -1.0;
  fTrueNCLepE = -1.0;
  fNetworkNCLepE = -1.0;

  fTrueNumHits = -1;
  fNetworkNumHits = -1;
}

DEFINE_ART_MODULE(extrapolation::RecoDump)
