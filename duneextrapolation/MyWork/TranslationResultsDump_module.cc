////////////////////////////////////////////////////////////////////////
// Class:       TranslationResultsDump
// Plugin Type: analyzer (Unknown Unknown)
// File:        TranslationResultsDump_module.cc
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
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

#include <string>
#include <map>

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace extrapolation {
  class TranslationResultsDump;
}


class extrapolation::TranslationResultsDump : public art::EDAnalyzer {
public:
  explicit TranslationResultsDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TranslationResultsDump(TranslationResultsDump const&) = delete;
  TranslationResultsDump(TranslationResultsDump&&) = delete;
  TranslationResultsDump& operator=(TranslationResultsDump const&) = delete;
  TranslationResultsDump& operator=(TranslationResultsDump&&) = delete;

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
  std::string fNDCVNResultsLabel;
  std::string fAltCVNResultsLabel;
  std::string fTrueHitsLabel; 
  std::string fNetworkHitsLabel;
  std::string fNDHitsLabel;
  std::string fAltHitsLabel;
  std::string fTrueNumuEResultsLabel;
  std::string fNetworkNumuEResultsLabel;
  std::string fAltNumuEResultsLabel;
  std::string fTrueNCEResultsLabel;
  std::string fNetworkNCEResultsLabel;
  std::string fAltNCEResultsLabel;

   // Options from fcl 
  bool fLegacy;
  
  TTree*             fTreeCVNResults;
  int                fRun;
  int                fSubRun;
  int                fEventNum;
  // Flavour scores
  float              fTrueNumuScore;
  float              fNetworkNumuScore;
  float              fNDNumuScore;
  float              fAltNumuScore;
  float              fTrueNueScore;
  float              fNetworkNueScore;
  float              fNDNueScore;
  float              fAltNueScore;
  float              fTrueNCScore;
  float              fNetworkNCScore;
  float              fNDNCScore;
  float              fAltNCScore;
  float              fTrueNutauScore;
  float              fNetworkNutauScore;
  float              fNDNutauScore;
  float              fAltNutauScore;
  // Pred number of daughter particle types
  int                fTrueNumPions;
  int                fNetworkNumPions;
  int                fNDNumPions;
  int                fAltNumPions;
  int                fTrueNumProtons;
  int                fNetworkNumProtons;
  int                fNDNumProtons;
  int                fAltNumProtons;
  int                fTrueNumPizeros;
  int                fNetworkNumPizeros;
  int                fNDNumPizeros;
  int                fAltNumPizeros;
  int                fTrueNumNeutrons;
  int                fNetworkNumNeutrons;
  int                fNDNumNeutrons;
  int                fAltNumNeutrons;
  // Antinuetrino score
  float              fTrueAntiNuScore;
  float              fNetworkAntiNuScore;
  float              fNDAntiNuScore;
  float              fAltAntiNuScore;
  // Interaction type score
  float              fTrueQEScore;
  float              fNetworkQEScore;
  float              fNDQEScore;
  float              fAltQEScore;
  float              fTrueResScore;
  float              fNetworkResScore;
  float              fNDResScore;
  float              fAltResScore;
  float              fTrueDISScore;
  float              fNetworkDISScore;
  float              fNDDISScore;
  float              fAltDISScore;
  float              fTrueOtherScore;
  float              fNetworkOtherScore;
  float              fNDOtherScore;
  float              fAltOtherScore;
  // Hit information
  float              fTrueHitIntegralSumZ;
  float              fNetworkHitIntegralSumZ;
  float              fNDHitIntegralSumZ;
  float              fAltHitIntegralSumZ;
  float              fTrueHitIntegralSumU;
  float              fNetworkHitIntegralSumU;
  float              fNDHitIntegralSumU;
  float              fAltHitIntegralSumU;
  float              fTrueHitIntegralSumV;
  float              fNetworkHitIntegralSumV;
  float              fNDHitIntegralSumV;
  float              fAltHitIntegralSumV;
  // True collection view hit information
  std::vector<float> fTrueHitTickWidthsZ;
  std::vector<int>   fTrueHitNumPerChZ;
  std::vector<float> fTrueHitIntegralsZ;
  // Nu energy reco information
  double             fTrueNumuNuE;
  double             fNetworkNumuNuE;
  double             fAltNumuNuE;
  double             fTrueNumuHadE;
  double             fNetworkNumuHadE;
  double             fAltNumuHadE;
  double             fTrueNumuLepE;
  double             fNetworkNumuLepE;
  double             fAltNumuLepE;
  double             fTrueNCNuE;
  double             fNetworkNCNuE;
  double             fAltNCNuE;
  double             fTrueNCHadE;
  double             fNetworkNCHadE;
  double             fAltNCHadE;
  double             fTrueNCLepE;
  double             fNetworkNCLepE;
  double             fAltNCLepE;
};


extrapolation::TranslationResultsDump::TranslationResultsDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTrueCVNResultsLabel      (p.get<std::string> ("TrueCVNResultsLabel")),
    fNetworkCVNResultsLabel   (p.get<std::string> ("NetworkCVNResultsLabel")),
    fNDCVNResultsLabel        (p.get<std::string> ("NDCVNResultsLabel")),
    fAltCVNResultsLabel       (p.get<std::string> ("AltCVNResultsLabel")),
    fTrueHitsLabel            (p.get<std::string> ("TrueHitsLabel")),    
    fNetworkHitsLabel         (p.get<std::string> ("NetworkHitsLabel")),      
    fNDHitsLabel              (p.get<std::string> ("NDHitsLabel")),
    fAltHitsLabel             (p.get<std::string> ("AltHitsLabel")),
    fTrueNumuEResultsLabel    (p.get<std::string> ("TrueNumuEResultsLabel")),
    fNetworkNumuEResultsLabel (p.get<std::string> ("NetworkNumuEResultsLabel")),
    fAltNumuEResultsLabel     (p.get<std::string> ("AltNumuEResultsLabel")),
    fTrueNCEResultsLabel      (p.get<std::string> ("TrueNCEResultsLabel")),
    fNetworkNCEResultsLabel   (p.get<std::string> ("NetworkNCEResultsLabel")),
    fAltNCEResultsLabel       (p.get<std::string> ("AltNCEResultsLabel")),
    fLegacy                   (p.get<bool> ("Legacy"))
{
  consumes<std::vector<cvn::Result>>(fTrueCVNResultsLabel);
  consumes<std::vector<cvn::Result>>(fNetworkCVNResultsLabel);
  consumes<std::vector<cvn::Result>>(fNDCVNResultsLabel);

  consumes<std::vector<recob::Hit>>(fTrueHitsLabel);
  consumes<std::vector<recob::Hit>>(fNetworkHitsLabel);
  consumes<std::vector<recob::Hit>>(fNDHitsLabel);

  if (!fLegacy) {
    consumes<dune::EnergyRecoOutput>(fTrueNumuEResultsLabel); 
    consumes<dune::EnergyRecoOutput>(fNetworkNumuEResultsLabel); 
    consumes<dune::EnergyRecoOutput>(fAltNumuEResultsLabel); 
    consumes<dune::EnergyRecoOutput>(fTrueNCEResultsLabel); 
    consumes<dune::EnergyRecoOutput>(fNetworkNCEResultsLabel); 
    consumes<dune::EnergyRecoOutput>(fAltNCEResultsLabel); 

    consumes<std::vector<cvn::Result>>(fAltCVNResultsLabel);

    consumes<std::vector<recob::Hit>>(fAltHitsLabel);
  }

  art::ServiceHandle<art::TFileService> tfs;

  fTreeCVNResults = tfs->make<TTree>("Results", "Results");
  fTreeCVNResults->Branch("Run", &fRun, "run/I");
  fTreeCVNResults->Branch("SubRun", &fSubRun, "subrun/I");
  fTreeCVNResults->Branch("EventNum", &fEventNum, "eventnum/I");
  // Flavour scores
  fTreeCVNResults->Branch("TrueNumuScore", &fTrueNumuScore, "truenumuscore/F");
  fTreeCVNResults->Branch("NetworkNumuScore", &fNetworkNumuScore, "networknumuscore/F");
  fTreeCVNResults->Branch("NDNumuScore", &fNDNumuScore, "ndnumuscore/F");
  fTreeCVNResults->Branch("AltNumuScore", &fAltNumuScore, "altnumuscore/F");
  fTreeCVNResults->Branch("TrueNueScore", &fTrueNueScore, "truenuescore/F");
  fTreeCVNResults->Branch("NetworkNueScore", &fNetworkNueScore, "networknuescore/F");
  fTreeCVNResults->Branch("NDNueScore", &fNDNueScore, "ndnuescore/F");  
  fTreeCVNResults->Branch("AltNueScore", &fAltNueScore, "altnuescore/F");  
  fTreeCVNResults->Branch("TrueNCScore", &fTrueNCScore, "truencscore/F");
  fTreeCVNResults->Branch("NetworkNCScore", &fNetworkNCScore, "networkncscore/F");
  fTreeCVNResults->Branch("NDNCScore", &fNDNCScore, "ndncscore/F");
  fTreeCVNResults->Branch("AltNCScore", &fAltNCScore, "altncscore/F");
  fTreeCVNResults->Branch("TrueNutauScore", &fTrueNutauScore, "truenutauscore/F");
  fTreeCVNResults->Branch("NetworkNutauScore", &fNetworkNutauScore, "networknutauscore/F");
  fTreeCVNResults->Branch("NDNutauScore", &fNDNutauScore, "ndnutauscore/F");
  fTreeCVNResults->Branch("AltNutauScore", &fAltNutauScore, "altnutauscore/F");
  // Pred number of daughter particles
  fTreeCVNResults->Branch("TrueNumPions", &fTrueNumPions, "truenumpions/I");
  fTreeCVNResults->Branch("NetworkNumPions", &fNetworkNumPions, "networknumpions/I");
  fTreeCVNResults->Branch("NDNumPions", &fNDNumPions, "ndnumpions/I");
  fTreeCVNResults->Branch("AltNumPions", &fAltNumPions, "altnumpions/I");
  fTreeCVNResults->Branch("TrueNumProtons", &fTrueNumProtons, "truenumprotons/I");
  fTreeCVNResults->Branch("NetworkNumProtons", &fNetworkNumProtons, "networknumprotons/I");
  fTreeCVNResults->Branch("NDNumProtons", &fNDNumProtons, "ndnumprotons/I");
  fTreeCVNResults->Branch("AltNumProtons", &fAltNumProtons, "altnumprotons/I");
  fTreeCVNResults->Branch("TrueNumPizeros", &fTrueNumPizeros, "truenumpizeros/I");
  fTreeCVNResults->Branch("NetworkNumPizeros", &fNetworkNumPizeros, "networknumpizeros/I");
  fTreeCVNResults->Branch("NDNumPizeros", &fNDNumPizeros, "ndnumpizeros/I");
  fTreeCVNResults->Branch("AltNumPizeros", &fAltNumPizeros, "altnumpizeros/I");
  fTreeCVNResults->Branch("TrueNumNeutrons", &fTrueNumNeutrons, "truenumneutrons/I");
  fTreeCVNResults->Branch("NetworkNumNeutrons", &fNetworkNumNeutrons, "networknumneutrons/I");
  fTreeCVNResults->Branch("NDNumNeutrons", &fNDNumNeutrons, "ndnumneutrons/I");
  fTreeCVNResults->Branch("AltNumNeutrons", &fAltNumNeutrons, "altnumneutrons/I");
  // Antineutrino score
  fTreeCVNResults->Branch("TrueAntiNuScore", &fTrueAntiNuScore, "trueantinuscore/F");
  fTreeCVNResults->Branch("NetworkAntiNuScore", &fNetworkAntiNuScore, "networkantinuscore/F");
  fTreeCVNResults->Branch("NDAntiNuScore", &fNDAntiNuScore, "ndantinuscore/F");
  fTreeCVNResults->Branch("AltAntiNuScore", &fAltAntiNuScore, "altantinuscore/F");
  // Interaction Type scores
  fTreeCVNResults->Branch("TrueQEScore", &fTrueQEScore, "trueqescore/F");
  fTreeCVNResults->Branch("NetworkQEScore", &fNetworkQEScore, "networkqescore/F");
  fTreeCVNResults->Branch("NDQEScore", &fNDQEScore, "ndqescore/F");
  fTreeCVNResults->Branch("AltQEScore", &fAltQEScore, "altqescore/F");
  fTreeCVNResults->Branch("TrueResScore", &fTrueResScore, "trueresscore/F");
  fTreeCVNResults->Branch("NetworkResScore", &fNetworkResScore, "networkresscore/F");
  fTreeCVNResults->Branch("NDResScore", &fNDResScore, "ndresscore/F");
  fTreeCVNResults->Branch("AltResScore", &fAltResScore, "altresscore/F");
  fTreeCVNResults->Branch("TrueDISScore", &fTrueDISScore, "truedisscore/F");
  fTreeCVNResults->Branch("NetworkDISScore", &fNetworkDISScore, "networkdisscore/F");
  fTreeCVNResults->Branch("NDDISScore", &fNDDISScore, "nddisscore/F");
  fTreeCVNResults->Branch("AltDISScore", &fAltDISScore, "altdisscore/F");
  fTreeCVNResults->Branch("TrueOtherScore", &fTrueOtherScore, "trueotherscore/F");
  fTreeCVNResults->Branch("NetworkOtherScore", &fNetworkOtherScore, "networkotherscore/F");
  fTreeCVNResults->Branch("NDOtherScore", &fNDQEScore, "ndotherscore/F"); 
  fTreeCVNResults->Branch("AltOtherScore", &fAltQEScore, "altotherscore/F"); 
  // Hit information
  fTreeCVNResults->Branch("TrueHitIntegralSumZ", &fTrueHitIntegralSumZ, "truehitintegralsumz/F");
  fTreeCVNResults->Branch("NetworkHitIntegralSumZ", &fNetworkHitIntegralSumZ, "networkhitintegralsumz/F");
  fTreeCVNResults->Branch("NDHitIntegralSumZ", &fNDHitIntegralSumZ, "ndhitintegralsumz/F");
  fTreeCVNResults->Branch("AltHitIntegralSumZ", &fAltHitIntegralSumZ, "althitintegralsumz/F");
  fTreeCVNResults->Branch("TrueHitIntegralSumU", &fTrueHitIntegralSumU, "truehitintegralsumu/F");
  fTreeCVNResults->Branch("NetworkHitIntegralSumU", &fNetworkHitIntegralSumU, "networkhitintegralsumu/F");
  fTreeCVNResults->Branch("NDHitIntegralSumU", &fNDHitIntegralSumU, "ndhitintegralsumu/F");
  fTreeCVNResults->Branch("AltHitIntegralSumU", &fAltHitIntegralSumU, "althitintegralsumu/F");
  fTreeCVNResults->Branch("TrueHitIntegralSumV", &fTrueHitIntegralSumV, "truehitintegralsumv/F");
  fTreeCVNResults->Branch("NetworkHitIntegralSumV", &fNetworkHitIntegralSumV, "networkhitintegralsumv/F");
  fTreeCVNResults->Branch("NDHitIntegralSumV", &fNDHitIntegralSumV, "ndhitintegralsumv/F");
  fTreeCVNResults->Branch("AltHitIntegralSumV", &fAltHitIntegralSumV, "althitintegralsumv/F");
  // True collection plane hit information
  fTreeCVNResults->Branch("TrueHitTickWidthsZ", &fTrueHitTickWidthsZ);
  fTreeCVNResults->Branch("TrueHitNumPerChZ", &fTrueHitNumPerChZ);
  fTreeCVNResults->Branch("TrueHitIntegralsZ", &fTrueHitIntegralsZ);
  // Nu energy reco information
  fTreeCVNResults->Branch("TrueNumuNuE", &fTrueNumuNuE, "truenumunue/D");
  fTreeCVNResults->Branch("NetworkNumuNuE", &fNetworkNumuNuE, "networknumunue/D");
  fTreeCVNResults->Branch("AltNumuNuE", &fAltNumuNuE, "altnumunue/D");
  fTreeCVNResults->Branch("TrueNumuHadE", &fTrueNumuHadE, "truenumuhade/D");
  fTreeCVNResults->Branch("NetworkNumuHadE", &fNetworkNumuHadE, "networknumuhade/D");
  fTreeCVNResults->Branch("AltNumuHadE", &fAltNumuHadE, "altnumuhade/D");
  fTreeCVNResults->Branch("TrueNumuLepE", &fTrueNumuLepE, "truenumulepe/D");
  fTreeCVNResults->Branch("NetworkNumuLepE", &fNetworkNumuLepE, "networknumulepe/D");
  fTreeCVNResults->Branch("AltNumuLepE", &fAltNumuLepE, "altnumulepe/D");
  fTreeCVNResults->Branch("TrueNCNuE", &fTrueNCNuE, "truencnue/D");
  fTreeCVNResults->Branch("NetworkNCNuE", &fNetworkNCNuE, "networkncnue/D");
  fTreeCVNResults->Branch("AltNCNuE", &fAltNCNuE, "altncnue/D");
  fTreeCVNResults->Branch("TrueNCHadE", &fTrueNCHadE, "truenchade/D");
  fTreeCVNResults->Branch("NetworkNCHadE", &fNetworkNCHadE, "networknchade/D");
  fTreeCVNResults->Branch("AltNCHadE", &fAltNCHadE, "altnchade/D");
  fTreeCVNResults->Branch("TrueNCLepE", &fTrueNCLepE, "truenclepe/D");
  fTreeCVNResults->Branch("NetworkNCLepE", &fNetworkNCLepE, "networknclepe/D");
  fTreeCVNResults->Branch("AltNCLepE", &fAltNCLepE, "altnclepe/D");
}

void extrapolation::TranslationResultsDump::analyze(art::Event const& e)
{
  this->reset();

  // Get event id info
  fRun = e.id().run();
  fSubRun = e.id().subRun(); 
  fEventNum = e.id().event();

  // Get results from CVN
  const auto trueCVNResults = e.getValidHandle<std::vector<cvn::Result>> (fTrueCVNResultsLabel);
  const auto networkCVNResults = e.getValidHandle<std::vector<cvn::Result>> (fNetworkCVNResultsLabel);
  const auto NDCVNResults = e.getValidHandle<std::vector<cvn::Result>> (fNDCVNResultsLabel);
  
  // Get hits
  const auto trueHits = e.getValidHandle<std::vector<recob::Hit>> (fTrueHitsLabel);
  const auto networkHits = e.getValidHandle<std::vector<recob::Hit>> (fNetworkHitsLabel);
  const auto NDHits = e.getValidHandle<std::vector<recob::Hit>> (fNDHitsLabel);
  
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
  
  // Get Hit information
  std::map<raw::ChannelID_t, int> chZHitCntr;
  for (const recob::Hit& hit : *trueHits) {
    if (hit.View() == geo::kZ) {
      fTrueHitIntegralSumZ += hit.Integral();
      fTrueHitTickWidthsZ.push_back(hit.EndTick() - hit.StartTick());
      fTrueHitIntegralsZ.push_back(hit.Integral());
      chZHitCntr[hit.Channel()]++;
    }
    else if (hit.View() == geo::kU) {
      fTrueHitIntegralSumU += hit.Integral();
    }
    else if (hit.View() == geo::kV) {
      fTrueHitIntegralSumV += hit.Integral();
    }
  }

  for (const recob::Hit& hit : *networkHits) {
    if (hit.View() == geo::kZ) {
      fNetworkHitIntegralSumZ += hit.Integral();
    }
    else if (hit.View() == geo::kU) {
      fNetworkHitIntegralSumU += hit.Integral();
    }
    else if (hit.View() == geo::kV) {
      fNetworkHitIntegralSumV += hit.Integral();
    }
  }

  for (const recob::Hit& hit : *NDHits) {
    if (hit.View() == geo::kZ) {
      fNDHitIntegralSumZ += hit.Integral();
    }
    else if (hit.View() == geo::kU) {
      fNDHitIntegralSumU += hit.Integral();
    }
    else if (hit.View() == geo::kV) {
      fNDHitIntegralSumV += hit.Integral();
    }
  }
  
  for (const auto& chCntrPair : chZHitCntr) {
    fTrueHitNumPerChZ.push_back(chCntrPair.second);
  }

  if (fLegacy) {
    return;
  }

  // Now assume energy reco output and second detsim process
  const auto altCVNResults = e.getValidHandle<std::vector<cvn::Result>> (fAltCVNResultsLabel);
  const auto altHits = e.getValidHandle<std::vector<recob::Hit>> (fAltHitsLabel);

  // Alt CVN info
  fAltNumuScore = altCVNResults->at(0).GetNumuProbability();
  fAltNueScore = altCVNResults->at(0).GetNueProbability();
  fAltNCScore = altCVNResults->at(0).GetNCProbability();
  fAltNutauScore = altCVNResults->at(0).GetNutauProbability();

  fAltNumPions = (int)altCVNResults->at(0).PredictedPions();
  fAltNumProtons = (int)altCVNResults->at(0).PredictedProtons();
  fAltNumPizeros = (int)altCVNResults->at(0).PredictedPizeros();
  fAltNumNeutrons = (int)altCVNResults->at(0).PredictedNeutrons();

  fAltAntiNuScore = altCVNResults->at(0).GetIsAntineutrinoProbability();

  fAltQEScore = altCVNResults->at(0).GetQEProbability();
  fAltResScore = altCVNResults->at(0).GetResProbability();
  fAltDISScore = altCVNResults->at(0).GetDISProbability();
  fAltOtherScore = altCVNResults->at(0).GetOtherProbability();

  // Alt Hit information
  for (const recob::Hit& hit : *altHits) {
    if (hit.View() == geo::kZ) {
      fAltHitIntegralSumZ += hit.Integral();
    }
    else if (hit.View() == geo::kU) {
      fAltHitIntegralSumU += hit.Integral();
    }
    else if (hit.View() == geo::kV) {
      fAltHitIntegralSumV += hit.Integral();
    }
  }

  // Get nu reco information
  const auto trueNumuEOut = e.getValidHandle<dune::EnergyRecoOutput> (fTrueNumuEResultsLabel);
  const auto networkNumuEOut = e.getValidHandle<dune::EnergyRecoOutput> (fNetworkNumuEResultsLabel);
  const auto altNumuEOut = e.getValidHandle<dune::EnergyRecoOutput> (fAltNumuEResultsLabel);
  const auto trueNCEOut = e.getValidHandle<dune::EnergyRecoOutput> (fTrueNCEResultsLabel);
  const auto networkNCEOut = e.getValidHandle<dune::EnergyRecoOutput> (fNetworkNCEResultsLabel);
  const auto altNCEOut = e.getValidHandle<dune::EnergyRecoOutput> (fAltNCEResultsLabel);

  fTrueNumuNuE = trueNumuEOut->fNuLorentzVector.E(); 
  fNetworkNumuNuE = networkNumuEOut->fNuLorentzVector.E(); 
  fAltNumuNuE = altNumuEOut->fNuLorentzVector.E(); 
  fTrueNumuHadE = trueNumuEOut->fHadLorentzVector.E(); 
  fNetworkNumuHadE = networkNumuEOut->fHadLorentzVector.E(); 
  fAltNumuHadE = altNumuEOut->fHadLorentzVector.E(); 
  fTrueNumuLepE = trueNumuEOut->fLepLorentzVector.E(); 
  fNetworkNumuLepE = networkNumuEOut->fLepLorentzVector.E(); 
  fAltNumuLepE = altNumuEOut->fLepLorentzVector.E(); 
  
  fTrueNCNuE = trueNCEOut->fNuLorentzVector.E(); 
  fNetworkNCNuE = networkNCEOut->fNuLorentzVector.E(); 
  fAltNCNuE = altNCEOut->fNuLorentzVector.E(); 
  fTrueNCHadE = trueNCEOut->fHadLorentzVector.E(); 
  fNetworkNCHadE = networkNCEOut->fHadLorentzVector.E(); 
  fAltNCHadE = altNCEOut->fHadLorentzVector.E(); 
  fTrueNCLepE = trueNCEOut->fLepLorentzVector.E(); 
  fNetworkNCLepE = networkNCEOut->fLepLorentzVector.E(); 
  fAltNCLepE = networkNCEOut->fLepLorentzVector.E(); 

  fTreeCVNResults->Fill();
}

void extrapolation::TranslationResultsDump::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
}

void extrapolation::TranslationResultsDump::endJob()
{
}

void extrapolation::TranslationResultsDump::reset() 
{
  fRun = -999;
  fSubRun = -999;
  fEventNum = -999;

  fTrueNumuScore = -999.0;
  fNetworkNumuScore = -999.0;
  fNDNumuScore = -999.0;
  fAltNumuScore = -999.0;
  fTrueNueScore = -999.0;
  fNetworkNueScore = -999.0;
  fNDNueScore = -999.0;
  fAltNueScore = -999.0;
  fTrueNCScore = -999.0;
  fNetworkNCScore = -999.0;
  fNDNCScore = -999.0;
  fAltNCScore = -999.0;
  fTrueNutauScore = -999.0;
  fNetworkNutauScore = -999.0;
  fNDNutauScore = -999.0;
  fAltNutauScore = -999.0;

  fTrueNumPions = -999;
  fNetworkNumPions = -999;
  fNDNumPions = -999;
  fAltNumPions = -999;
  fTrueNumProtons = -999;
  fNetworkNumProtons = -999;
  fNDNumProtons = -999;
  fAltNumProtons = -999;
  fTrueNumPizeros = -999;
  fNetworkNumPizeros = -999;
  fNDNumPizeros = -999;
  fAltNumPizeros = -999;
  fTrueNumNeutrons = -999;
  fNetworkNumNeutrons = -999;
  fNDNumNeutrons = -999;
  fAltNumNeutrons = -999;

  fTrueAntiNuScore = -999.0;
  fNetworkAntiNuScore = -999.0;
  fNDAntiNuScore = -999.0;
  fAltAntiNuScore = -999.0;

  fTrueQEScore = -999.0;
  fNetworkQEScore = -999.0;
  fNDQEScore = -999.0;
  fAltQEScore = -999.0;
  fTrueResScore = -999.0;
  fNetworkResScore = -999.0;
  fNDResScore = -999.0;
  fAltResScore = -999.0;
  fTrueDISScore = -999.0;
  fNetworkDISScore = -999.0;
  fNDDISScore = -999.0;
  fAltDISScore = -999.0;
  fTrueOtherScore = -999.0;
  fNetworkOtherScore = -999.0;
  fNDOtherScore = -999.0;
  fAltOtherScore = -999.0;
  
  fTrueHitIntegralSumZ = 0.0;
  fNetworkHitIntegralSumZ = 0.0;
  fNDHitIntegralSumZ = 0.0;
  fAltHitIntegralSumZ = 0.0;
  fTrueHitIntegralSumU = 0.0;
  fNetworkHitIntegralSumU = 0.0;
  fNDHitIntegralSumU = 0.0;
  fAltHitIntegralSumU = 0.0;
  fTrueHitIntegralSumV = 0.0;
  fNetworkHitIntegralSumV = 0.0;
  fNDHitIntegralSumV = 0.0;
  fAltHitIntegralSumV = 0.0;

  fTrueHitTickWidthsZ.clear();
  fTrueHitNumPerChZ.clear();
  fTrueHitIntegralsZ.clear();

  fTrueNumuNuE = -1.0;
  fNetworkNumuNuE = -1.0;
  fAltNumuNuE = -1.0;
  fTrueNumuHadE = -1.0;
  fNetworkNumuHadE = -1.0;
  fAltNumuHadE = -1.0;
  fTrueNumuLepE = -1.0;
  fNetworkNumuLepE = -1.0;
  fAltNumuLepE = -1.0;
  fTrueNCNuE = -1.0;
  fNetworkNCNuE = -1.0;
  fAltNCNuE = -1.0;
  fTrueNCHadE = -1.0;
  fNetworkNCHadE = -1.0;
  fAltNCHadE = -1.0;
  fTrueNCLepE = -1.0;
  fNetworkNCLepE = -1.0;
  fAltNCLepE = -1.0;
}

DEFINE_ART_MODULE(extrapolation::TranslationResultsDump)
