////////////////////////////////////////////////////////////////////////
// Class:       AddFDReco
// Plugin Type: analyzer (Unknown Unknown)
// File:        AddFDReco.cc
//
// Crated on 20 Sep 23 Alex Wilkinson
// Dump reco data products to the larnd-sim HDF5 file
// Made for creating a dataset of paired ND det resp and FD reco
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
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "highfive/H5DataSet.hpp"
#include "highfive/H5File.hpp"
#include "highfive/H5Group.hpp"
#include "highfive/H5Object.hpp"
#include "highfive/H5DataType.hpp"
#include "highfive/H5DataSpace.hpp"

#include <string>
#include <vector>
#include <map>

typedef struct recoFD {
  int eventID;
  float numu_score;
  float nue_score;
  float nuc_score;
  float nutau_score;
  float antinu_score;
  float p_0_score;
  float p_1_score;
  float p_2_score;
  float p_N_score;
  float pi_0_score;
  float pi_1_score;
  float pi_2_score;
  float pi_N_score;
  float pi0_0_score;
  float pi0_1_score;
  float pi0_2_score;
  float pi0_N_score;
  float n_0_score;
  float n_1_score;
  float n_2_score;
  float n_N_score;
  float numu_nu_E;
  float numu_had_E;
  float numu_lep_E;
  int numu_reco_method;
  int numu_longest_track_contained;
  int numu_longest_track_mom_method;
  float nue_nu_E;
  float nue_had_E;
  float nue_lep_E;
  int nue_reco_method;
  float nc_nu_E;
  float nc_had_E;
  float nc_lep_E;
  int nc_reco_method;
} recoFD;

HighFive::CompoundType make_recoFD() {
  return {
    {"eventID", HighFive::AtomicType<int>{}},
    {"numu_score", HighFive::AtomicType<float>{}},
    {"nue_score", HighFive::AtomicType<float>{}},
    {"nuc_score", HighFive::AtomicType<float>{}},
    {"nutau_score", HighFive::AtomicType<float>{}},
    {"antinu_score", HighFive::AtomicType<float>{}},
    {"0_p_score", HighFive::AtomicType<float>{}},
    {"1_p_score", HighFive::AtomicType<float>{}},
    {"2_p_score", HighFive::AtomicType<float>{}},
    {"N_p_score", HighFive::AtomicType<float>{}},
    {"0_pi_score", HighFive::AtomicType<float>{}},
    {"1_pi_score", HighFive::AtomicType<float>{}},
    {"2_pi_score", HighFive::AtomicType<float>{}},
    {"N_pi_score", HighFive::AtomicType<float>{}},
    {"0_pi0_score", HighFive::AtomicType<float>{}},
    {"1_pi0_score", HighFive::AtomicType<float>{}},
    {"2_pi0_score", HighFive::AtomicType<float>{}},
    {"N_pi0_score", HighFive::AtomicType<float>{}},
    {"0_n_score", HighFive::AtomicType<float>{}},
    {"1_n_score", HighFive::AtomicType<float>{}},
    {"2_n_score", HighFive::AtomicType<float>{}},
    {"N_n_score", HighFive::AtomicType<float>{}},
    {"numu_nu_E", HighFive::AtomicType<float>{}},
    {"numu_had_E", HighFive::AtomicType<float>{}},
    {"numu_lep_E", HighFive::AtomicType<float>{}},
    {"numu_reco_method", HighFive::AtomicType<int>{}},
    {"numu_longest_track_contained", HighFive::AtomicType<int>{}},
    {"numu_longest_track_mom_method", HighFive::AtomicType<int>{}},
    {"nue_nu_E", HighFive::AtomicType<float>{}},
    {"nue_had_E", HighFive::AtomicType<float>{}},
    {"nue_lep_E", HighFive::AtomicType<float>{}},
    {"nue_reco_method", HighFive::AtomicType<int>{}},
    {"nc_nu_E", HighFive::AtomicType<float>{}},
    {"nc_had_E", HighFive::AtomicType<float>{}},
    {"nc_lep_E", HighFive::AtomicType<float>{}},
    {"nc_reco_method", HighFive::AtomicType<int>{}}
  };
}

HIGHFIVE_REGISTER_TYPE(recoFD, make_recoFD)

namespace extrapolation {
  class AddFDReco;
}

class extrapolation::AddFDReco : public art::EDAnalyzer {
public:
  explicit AddFDReco(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  AddFDReco(AddFDReco const&) = delete;
  AddFDReco(AddFDReco&&) = delete;
  AddFDReco& operator=(AddFDReco const&) = delete;
  AddFDReco& operator=(AddFDReco&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Methods
  recoFD getReco(
    art::Event const& e,
    int eventID,
    std::string& CVNResultsLabel,
    std::string& numuEResultsLabel,
    std::string& nueEResultsLabel,
    std::string& NCEResultsLabel
  );

  // Members
  const geo::GeometryCore* fGeom;

  HighFive::File* fFile;

  // Data to write out to hdf5
  std::vector<recoFD> fReco;

  // Product labels
  std::string fEventIDSEDLabel;
  std::string fCVNResultsLabel;
  std::string fNumuEResultsLabel;
  std::string fNueEResultsLabel;
  std::string fNCEResultsLabel;

  std::string fNDFDH5FileLoc;
};


extrapolation::AddFDReco::AddFDReco(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fEventIDSEDLabel   (p.get<std::string>("EventIDSEDLabel")),
    fCVNResultsLabel   (p.get<std::string>("CVNResultsLabel")),
    fNumuEResultsLabel (p.get<std::string>("NumuEResultsLabel")),
    fNueEResultsLabel  (p.get<std::string>("NueEResultsLabel")),
    fNCEResultsLabel   (p.get<std::string>("NCEResultsLabel")),
    fNDFDH5FileLoc       (p.get<std::string>("NDFDH5FileLoc"))
{
  consumes<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);

  consumes<std::vector<cvn::Result>>(fCVNResultsLabel);

  consumes<dune::EnergyRecoOutput>(fNumuEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNueEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNCEResultsLabel);
}

void extrapolation::AddFDReco::analyze(art::Event const& e)
{
  // Get eventID
  const auto eventIDSED = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);
  int eventID = (*eventIDSED)[0].TrackID();

  // Store reco for this event
  recoFD eventReco = getReco(
    e, eventID, fCVNResultsLabel, fNumuEResultsLabel, fNueEResultsLabel, fNCEResultsLabel
  );
  fReco.push_back(eventReco);
}

void extrapolation::AddFDReco::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  fFile = new HighFive::File(fNDFDH5FileLoc, HighFive::File::ReadWrite);
}

void extrapolation::AddFDReco::endJob()
{
  // Write reco to hdf5
  fFile->createDataSet("fd_reco", fReco);
}

recoFD extrapolation::AddFDReco::getReco(
  art::Event const& e,
  int eventID,
  std::string& CVNResultsLabel,
  std::string& numuEResultsLabel,
  std::string& nueEResultsLabel,
  std::string& NCEResultsLabel
)
{
  // Get results from CVN
  const auto CVNResults = e.getValidHandle<std::vector<cvn::Result>>(CVNResultsLabel);
  // Happens rarely, I guess its when the pixel map producer fails (too few hits maybe?)
  const bool validCVN = CVNResults->size() != 0;

  // Get flavour scores
  float numuScore = validCVN ? CVNResults->at(0).GetNumuProbability() : -1;
  float nueScore = validCVN ? CVNResults->at(0).GetNueProbability() : -1;
  float ncScore = validCVN ? CVNResults->at(0).GetNCProbability() : -1;
  float nutauScore = validCVN ? CVNResults->at(0).GetNutauProbability() : -1;

  // Get antineutrino score
  float antiNuScore = validCVN ? CVNResults->at(0).GetIsAntineutrinoProbability() : -1;

  // Get CVN doodads
  float proton0Score = validCVN ? CVNResults->at(0).Get0protonsProbability() : -1;
  float proton1Score = validCVN ? CVNResults->at(0).Get1protonsProbability() : -1;
  float proton2Score = validCVN ? CVNResults->at(0).Get2protonsProbability() : -1;
  float protonNScore = validCVN ? CVNResults->at(0).GetNprotonsProbability() : -1;
  float pion0Score = validCVN ? CVNResults->at(0).Get0pionsProbability() : -1;
  float pion1Score = validCVN ? CVNResults->at(0).Get1pionsProbability() : -1;
  float pion2Score = validCVN ? CVNResults->at(0).Get2pionsProbability() : -1;
  float pionNScore = validCVN ? CVNResults->at(0).GetNpionsProbability() : -1;
  float pionZero0Score = validCVN ? CVNResults->at(0).Get0pizerosProbability() : -1;
  float pionZero1Score = validCVN ? CVNResults->at(0).Get1pizerosProbability() : -1;
  float pionZero2Score = validCVN ? CVNResults->at(0).Get2pizerosProbability() : -1;
  float pionZeroNScore = validCVN ? CVNResults->at(0).GetNpizerosProbability() : -1;
  float neutron0Score = validCVN ? CVNResults->at(0).Get0neutronsProbability() : -1;
  float neutron1Score = validCVN ? CVNResults->at(0).Get1neutronsProbability() : -1;
  float neutron2Score = validCVN ? CVNResults->at(0).Get2neutronsProbability() : -1;
  float neutronNScore = validCVN ? CVNResults->at(0).GetNneutronsProbability() : -1;

  // Get nu reco information
  const auto numuEOut = e.getValidHandle<dune::EnergyRecoOutput>(numuEResultsLabel);
  const auto nueEOut = e.getValidHandle<dune::EnergyRecoOutput>(nueEResultsLabel);
  const auto NCEOut = e.getValidHandle<dune::EnergyRecoOutput>(NCEResultsLabel);

  float numuNuE = (float)numuEOut->fNuLorentzVector.E();
  float numuHadE = (float)numuEOut->fHadLorentzVector.E();
  float numuLepE = (float)numuEOut->fLepLorentzVector.E();
  int numuRecoMethod = (int)numuEOut->recoMethodUsed;
  int numuLongestTrackContained = (int)numuEOut->longestTrackContained;
  int numuTrackMomMethod = (int)numuEOut->trackMomMethod;

  float nueNuE = (float)nueEOut->fNuLorentzVector.E();
  float nueHadE = (float)nueEOut->fHadLorentzVector.E();
  float nueLepE = (float)nueEOut->fLepLorentzVector.E();
  int nueRecoMethod = (int)nueEOut->recoMethodUsed;

  float ncNuE = (float)NCEOut->fNuLorentzVector.E();
  float ncHadE = (float)NCEOut->fHadLorentzVector.E();
  float ncLepE = (float)NCEOut->fLepLorentzVector.E();
  int ncRecoMethod = (int)NCEOut->recoMethodUsed;

  recoFD eventReco = {
    eventID,
    numuScore,
    nueScore,
    ncScore,
    nutauScore,
    antiNuScore,
    proton0Score,
    proton1Score,
    proton2Score,
    protonNScore,
    pion0Score,
    pion1Score,
    pion2Score,
    pionNScore,
    pionZero0Score,
    pionZero1Score,
    pionZero2Score,
    pionZeroNScore,
    neutron0Score,
    neutron1Score,
    neutron2Score,
    neutronNScore,
    numuNuE,
    numuHadE,
    numuLepE,
    numuRecoMethod,
    numuLongestTrackContained,
    numuTrackMomMethod,
    nueNuE,
    nueHadE,
    nueLepE,
    nueRecoMethod,
    ncNuE,
    ncHadE,
    ncLepE,
    ncRecoMethod
  };

  return eventReco;
}

DEFINE_ART_MODULE(extrapolation::AddFDReco)
