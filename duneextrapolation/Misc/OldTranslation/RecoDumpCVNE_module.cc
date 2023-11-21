////////////////////////////////////////////////////////////////////////
// Class:       RecoDumpCVNECVNE
// Plugin Type: analyzer (Unknown Unknown)
// File:        RecoDumpCVNECVNE_module.cc
//
// Crated on 12 May 23 Alex Wilkinson
// Dump reco data products to a flat TTree
// Made for creating a dataset of paired ND and FD reconstruction
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
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include <string>
#include <vector>
#include <map>

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace extrapolation {
  class RecoDumpCVNE;
}


class extrapolation::RecoDumpCVNE : public art::EDAnalyzer {
public:
  explicit RecoDumpCVNE(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RecoDumpCVNE(RecoDumpCVNE const&) = delete;
  RecoDumpCVNE(RecoDumpCVNE&&) = delete;
  RecoDumpCVNE& operator=(RecoDumpCVNE const&) = delete;
  RecoDumpCVNE& operator=(RecoDumpCVNE&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My function.
  void reset();

private:
  const geo::GeometryCore* fGeom;

  bool                 fCalcOldHadEReco;
  calo::CalorimetryAlg fCaloAlg;
  double               fRecombFactor;
  // Labels from fcl
  std::string fCVNResultsLabel;
  std::string fNumuEResultsLabel;
  std::string fNueEResultsLabel;
  std::string fNCEResultsLabel;
  std::string fEventIDSEDLabel;
  std::string fOldHadERecoHitLabel;

  TTree* fTreeReco;
  int    fEventID;
  // Flavour scores
  float fNumuScore;
  float fNueScore;
  float fNCScore;
  float fNutauScore;
  // Antinuetrino score
  float fAntiNuScore;
  // Other CVN doodads
  float f0ProtonScore;
  float f1ProtonScore;
  float f2ProtonScore;
  float fNProtonScore;
  float f0PionScore;
  float f1PionScore;
  float f2PionScore;
  float fNPionScore;
  float f0PionZeroScore;
  float f1PionZeroScore;
  float f2PionZeroScore;
  float fNPionZeroScore;
  float f0NeutronScore;
  float f1NeutronScore;
  float f2NeutronScore;
  float fNNeutronScore;
  // Nu energy reco information
  float fNumuNuE;
  float fNumuLepE;
  float fNumuHadE;
  int   fNumuRecoMethod;
  int   fNumuLongestTrackContained;
  int   fNumuTrackMomMethod;
  float fNueNuE;
  float fNueLepE;
  float fNueHadE;
  int   fNueRecoMethod;
  float fNCNuE;
  float fNCLepE;
  float fNCHadE;
  int   fNCRecoMethod;
  // Nu Had energy reco old implementation
  float fEDepP;
  float fEDepN;
  float fEDepPip;
  float fEDepPim;
  float fEDepPi0;
  float fEDepOther;

  // Protect against indexing empty reco vectors
  int fNumOORErrs;
};


extrapolation::RecoDumpCVNE::RecoDumpCVNE(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fCalcOldHadEReco     (p.get<bool>("CalcOldHadEReco")),
    fCaloAlg             (p.get<fhicl::ParameterSet>("CaloAlg")),
    fRecombFactor        (p.get<double>("RecombFactor")),
    fCVNResultsLabel     (p.get<std::string>("CVNResultsLabel")),
    fNumuEResultsLabel   (p.get<std::string>("NumuEResultsLabel")),
    fNueEResultsLabel    (p.get<std::string>("NueEResultsLabel")),
    fNCEResultsLabel     (p.get<std::string>("NCEResultsLabel")),
    fEventIDSEDLabel     (p.get<std::string>("EventIDSEDLabel")),
    fOldHadERecoHitLabel (p.get<std::string>("OldHadRecoHitLabel"))

{
  consumes<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);

  consumes<std::vector<cvn::Result>>(fCVNResultsLabel);

  consumes<dune::EnergyRecoOutput>(fNumuEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNueEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNCEResultsLabel);

  consumes<std::vector<recob::Hit>>(fNCEResultsLabel);

  art::ServiceHandle<art::TFileService> tfs;

  fTreeReco = tfs->make<TTree>("FDReco", "FDReco");
  fTreeReco->Branch("eventID", &fEventID, "eventid/I");
  // Flavour scores
  fTreeReco->Branch("numuScore", &fNumuScore, "numuScore/F");
  fTreeReco->Branch("nueScore", &fNueScore, "nueScore/F");
  fTreeReco->Branch("NCScore", &fNCScore, "NCScore/F");
  fTreeReco->Branch("nutauScore", &fNutauScore, "nutauScore/F");
  // Antineutrino score
  fTreeReco->Branch("antiNuScore", &fAntiNuScore, "antiNuScore/F");
  // Other CVN doodada
  fTreeReco->Branch("proton0Score", &f0ProtonScore, "proton0Score/F");
  fTreeReco->Branch("proton1Score", &f1ProtonScore, "proton1Score/F");
  fTreeReco->Branch("proton2Score", &f2ProtonScore, "proton2Score/F");
  fTreeReco->Branch("protonNScore", &fNProtonScore, "protonNScore/F");
  fTreeReco->Branch("pion0Score", &f0PionScore, "pion0Score/F");
  fTreeReco->Branch("pion1Score", &f1PionScore, "pion1Score/F");
  fTreeReco->Branch("pion2Score", &f2PionScore, "pion2Score/F");
  fTreeReco->Branch("pionNScore", &fNPionScore, "pionNScore/F");
  fTreeReco->Branch("pionZero0Score", &f0PionZeroScore, "pionZero0Score/F");
  fTreeReco->Branch("pionZero1Score", &f1PionZeroScore, "pionZero1Score/F");
  fTreeReco->Branch("pionZero2Score", &f2PionZeroScore, "pionZero2Score/F");
  fTreeReco->Branch("pionZeroNScore", &fNPionZeroScore, "pionZeroNScore/F");
  fTreeReco->Branch("neutron0Score", &f0NeutronScore, "neutron0Score/F");
  fTreeReco->Branch("neutron1Score", &f1NeutronScore, "neutron1Score/F");
  fTreeReco->Branch("neutron2Score", &f2NeutronScore, "neutron2Score/F");
  fTreeReco->Branch("neutronNScore", &fNNeutronScore, "neutronNScore/F");
  // Nu energy reco information
  fTreeReco->Branch("numuNuE", &fNumuNuE, "numuNuE/F");
  fTreeReco->Branch("numuHadE", &fNumuHadE, "numuHadE/F");
  fTreeReco->Branch("numuLepE", &fNumuLepE, "numuLepE/F");
  fTreeReco->Branch("numuRecoMethod", &fNumuRecoMethod, "numuRecoMethod/I");
  fTreeReco->Branch(
    "numuLongestTrackContained", &fNumuLongestTrackContained, "numuLongestTrackContained/I");
  fTreeReco->Branch("numuTrackMomMethod", &fNumuTrackMomMethod, "numuTrackMomMethod/I");
  fTreeReco->Branch("nueNuE", &fNueNuE, "nueNuE/F");
  fTreeReco->Branch("nueHadE", &fNueHadE, "nueHadE/F");
  fTreeReco->Branch("nueLepE", &fNueLepE, "nueLepE/F");
  fTreeReco->Branch("nueRecoMethod", &fNueRecoMethod, "nueRecoMethod/I");
  fTreeReco->Branch("NCNuE", &fNCNuE, "NCNuE/F");
  fTreeReco->Branch("NCHadE", &fNCHadE, "NCHadE/F");
  fTreeReco->Branch("NCLepE", &fNCLepE, "NCLepE/F");
  fTreeReco->Branch("NCRecoMethod", &fNCRecoMethod, "NCRecoMethod/I");
  // Nu Had energy reco old implementation
  if (fCalcOldHadEReco) {
    fTreeReco->Branch("eDepP", &fEDepP, "eDepP/F");
    fTreeReco->Branch("eDepN", &fEDepN, "eDepN/F");
    fTreeReco->Branch("eDepPip", &fEDepPip, "eDepPip/F");
    fTreeReco->Branch("eDepPim", &fEDepPim, "eDepPim/F");
    fTreeReco->Branch("eDepPi0", &fEDepPi0, "eDepPi0/F");
    fTreeReco->Branch("eDepOther", &fEDepOther, "eDepOther/F");
  }
}

void extrapolation::RecoDumpCVNE::analyze(art::Event const& e)
{
  this->reset();

  try {
    // Get SED containing eventID
    const auto eventIDSEDs =
      e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);
    fEventID = eventIDSEDs->at(0).TrackID();

    // Get results from CVN
    const auto CVNResults = e.getValidHandle<std::vector<cvn::Result>>(fCVNResultsLabel);

    // Get flavour scores
    fNumuScore = CVNResults->at(0).GetNumuProbability();
    fNueScore = CVNResults->at(0).GetNueProbability();
    fNCScore = CVNResults->at(0).GetNCProbability();
    fNutauScore = CVNResults->at(0).GetNutauProbability();

    // Get antineutrino score
    fAntiNuScore = CVNResults->at(0).GetIsAntineutrinoProbability();

    // Get CVN doodads
    f0ProtonScore = CVNResults->at(0).Get0protonsProbability();
    f1ProtonScore = CVNResults->at(0).Get1protonsProbability();
    f2ProtonScore = CVNResults->at(0).Get2protonsProbability();
    fNProtonScore = CVNResults->at(0).GetNprotonsProbability();
    f0PionScore = CVNResults->at(0).Get0pionsProbability();
    f1PionScore = CVNResults->at(0).Get1pionsProbability();
    f2PionScore = CVNResults->at(0).Get2pionsProbability();
    fNPionScore = CVNResults->at(0).GetNpionsProbability();
    f0PionZeroScore = CVNResults->at(0).Get0pizerosProbability();
    f1PionZeroScore = CVNResults->at(0).Get1pizerosProbability();
    f2PionZeroScore = CVNResults->at(0).Get2pizerosProbability();
    fNPionZeroScore = CVNResults->at(0).GetNpizerosProbability();
    f0NeutronScore = CVNResults->at(0).Get0neutronsProbability();
    f1NeutronScore = CVNResults->at(0).Get1neutronsProbability();
    f2NeutronScore = CVNResults->at(0).Get2neutronsProbability();
    fNNeutronScore = CVNResults->at(0).GetNneutronsProbability();

    // Get nu reco information
    const auto numuEOut = e.getValidHandle<dune::EnergyRecoOutput>(fNumuEResultsLabel);
    const auto nueEOut = e.getValidHandle<dune::EnergyRecoOutput>(fNueEResultsLabel);
    const auto NCEOut = e.getValidHandle<dune::EnergyRecoOutput>(fNCEResultsLabel);

    fNumuNuE = (float)numuEOut->fNuLorentzVector.E();
    fNumuHadE = (float)numuEOut->fHadLorentzVector.E();
    fNumuLepE = (float)numuEOut->fLepLorentzVector.E();
    fNumuRecoMethod = (int)numuEOut->recoMethodUsed;
    fNumuLongestTrackContained = (int)numuEOut->longestTrackContained;
    fNumuTrackMomMethod = (int)numuEOut->trackMomMethod;

    fNueNuE = (float)nueEOut->fNuLorentzVector.E();
    fNueHadE = (float)nueEOut->fHadLorentzVector.E();
    fNueLepE = (float)nueEOut->fLepLorentzVector.E();
    fNueRecoMethod = (int)nueEOut->recoMethodUsed;

    fNCNuE = (float)NCEOut->fNuLorentzVector.E();
    fNCHadE = (float)NCEOut->fHadLorentzVector.E();
    fNCLepE = (float)NCEOut->fLepLorentzVector.E();
    fNCRecoMethod = (int)NCEOut->recoMethodUsed;

    if (fCalcOldHadEReco) {
      // auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
      // double t0 = detprop->TriggerOffset();

      // std::map<int, double> tIDEDep;

      // const auto hits = e.getValidHandle<std::vector<recob::Hit>>(fOldHadERecoHitLabel);

      // for (recob::Hit hit : *hits) {
      //   double qeLifetimeCorrected = hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
      // }

      // NOTE stopped implementing as I think there is no point trying to replicate this visible
      // energy variable identically as FD sim/reco has changed so much since TDR days so that
      // would need to be addressed also
    }

    fTreeReco->Fill();
  }
  catch (const std::out_of_range& err) {
    fNumOORErrs++;

    this->reset();
    fTreeReco->Fill();
  }
}

void extrapolation::RecoDumpCVNE::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  fNumOORErrs = 0;
}

void extrapolation::RecoDumpCVNE::endJob()
{
  std::cout << "There were " << fNumOORErrs <<
               " out of range errors when accessing data products\n";
}

void extrapolation::RecoDumpCVNE::reset()
{
  fEventID = -99;

  fNumuScore = -99.0;
  fNueScore = -99.0;
  fNCScore = -99.0;
  fNutauScore = -99.0;
  fAntiNuScore = -99.0;

  f0ProtonScore = -99.0;
  f1ProtonScore = -99.0;
  f2ProtonScore = -99.0;
  fNProtonScore = -99.0;
  f0PionScore = -99.0;
  f1PionScore = -99.0;
  f2PionScore = -99.0;
  fNPionScore = -99.0;
  f0PionZeroScore = -99.0;
  f1PionZeroScore = -99.0;
  f2PionZeroScore = -99.0;
  fNPionZeroScore = -99.0;
  f0NeutronScore = -99.0;
  f1NeutronScore = -99.0;
  f2NeutronScore = -99.0;
  fNNeutronScore = -99.0;

  fNumuNuE = -99.0;
  fNumuLepE = -99.0;
  fNumuHadE = -99.0;
  fNumuRecoMethod = -99;
  fNumuLongestTrackContained = -99;
  fNumuTrackMomMethod = -99;
  fNueNuE = -99.0;
  fNueLepE = -99.0;
  fNueHadE = -99.0;
  fNueRecoMethod = -99;
  fNCNuE = -99.0;
  fNCLepE = -99.0;
  fNCHadE = -99.0;
  fNCRecoMethod = -99;

  fEDepP = -99.0;
  fEDepN = -99.0;
  fEDepPip = -99.0;
  fEDepPim = -99.0;
  fEDepPi0 = -99.0;
  fEDepOther = -99.0;
}

DEFINE_ART_MODULE(extrapolation::RecoDumpCVNE)
