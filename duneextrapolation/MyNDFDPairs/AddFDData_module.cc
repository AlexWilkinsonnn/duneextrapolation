////////////////////////////////////////////////////////////////////////
// Class:       AddFDData
// Plugin Type: analyzer (Unknown Unknown)
// File:        AddFDData_module.cc
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

typedef struct packet3d {
  int eventID;
  double adc;
  double x;
  double x_module;
  double y;
  double z;
  double z_module;
} packet3d;

HighFive::CompoundType make_packet3d() {
  return {
    {"eventID", HighFive::AtomicType<int>{}},
    {"adc", HighFive::AtomicType<double>{}},
    {"x", HighFive::AtomicType<double>{}},
    {"x_module", HighFive::AtomicType<double>{}},
    {"y", HighFive::AtomicType<double>{}},
    {"z", HighFive::AtomicType<double>{}},
    {"z_module", HighFive::AtomicType<double>{}}
  };
}

typedef struct packetProj {
  int eventID;
  float adc;
  int local_ch;
  int tick;
  float drift_dist;
} packetProj;

HighFive::CompoundType make_packetProj() {
  return {
    {"eventID", HighFive::AtomicType<int>{}},
    {"adc", HighFive::AtomicType<float>{}},
    {"local_ch", HighFive::AtomicType<int>{}},
    {"tick", HighFive::AtomicType<int>{}},
    {"drift_dist", HighFive::AtomicType<float>{}}
  };
}

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

HIGHFIVE_REGISTER_TYPE(packet3d, make_packet3d)
HIGHFIVE_REGISTER_TYPE(packetProj, make_packetProj)
HIGHFIVE_REGISTER_TYPE(recoFD, make_recoFD)

namespace extrapolation {
  class AddFDData;
}

class extrapolation::AddFDData : public art::EDAnalyzer {
public:
  explicit AddFDData(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  AddFDData(AddFDData const&) = delete;
  AddFDData(AddFDData&&) = delete;
  AddFDData& operator=(AddFDData const&) = delete;
  AddFDData& operator=(AddFDData&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My function.
  void reset();

private:
  const geo::GeometryCore* fGeom;

  HighFive::File* fFile;

  // std::vector<packetProj> fProjsZ;
  std::vector<recoFD> fReco;

  // Labels from fcl
  std::string fEventIDSEDLabel;
  std::string fCVNResultsLabel;
  std::string fNumuEResultsLabel;
  std::string fNueEResultsLabel;
  std::string fNCEResultsLabel;

  std::string fNDH5FileLoc;
};


extrapolation::AddFDData::AddFDData(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fEventIDSEDLabel   (p.get<std::string>("EventIDSEDLabel")),
    fCVNResultsLabel   (p.get<std::string>("CVNResultsLabel")),
    fNumuEResultsLabel (p.get<std::string>("NumuEResultsLabel")),
    fNueEResultsLabel  (p.get<std::string>("NueEResultsLabel")),
    fNCEResultsLabel   (p.get<std::string>("NCEResultsLabel")),
    fNDH5FileLoc       (p.get<std::string>("NDH5FileLoc"))
{
  consumes<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);

  consumes<std::vector<cvn::Result>>(fCVNResultsLabel);

  consumes<std::vector<cvn::Result>>(fCVNResultsLabel);

  consumes<dune::EnergyRecoOutput>(fNumuEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNueEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNCEResultsLabel);
}

void extrapolation::AddFDData::analyze(art::Event const& e)
{
  // Get eventID
  const auto eventIDSED = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);
  int eventID = (*eventIDSED)[0].TrackID();

  // Get results from CVN
  const auto CVNResults = e.getValidHandle<std::vector<cvn::Result>>(fCVNResultsLabel);

  // Get flavour scores
  float numuScore = CVNResults->at(0).GetNumuProbability();
  float nueScore = CVNResults->at(0).GetNueProbability();
  float ncScore = CVNResults->at(0).GetNCProbability();
  float nutauScore = CVNResults->at(0).GetNutauProbability();

  // Get antineutrino score
  float antiNuScore = CVNResults->at(0).GetIsAntineutrinoProbability();

  // Get CVN doodads
  float proton0Score = CVNResults->at(0).Get0protonsProbability();
  float proton1Score = CVNResults->at(0).Get1protonsProbability();
  float proton2Score = CVNResults->at(0).Get2protonsProbability();
  float protonNScore = CVNResults->at(0).GetNprotonsProbability();
  float pion0Score = CVNResults->at(0).Get0pionsProbability();
  float pion1Score = CVNResults->at(0).Get1pionsProbability();
  float pion2Score = CVNResults->at(0).Get2pionsProbability();
  float pionNScore = CVNResults->at(0).GetNpionsProbability();
  float pionZero0Score = CVNResults->at(0).Get0pizerosProbability();
  float pionZero1Score = CVNResults->at(0).Get1pizerosProbability();
  float pionZero2Score = CVNResults->at(0).Get2pizerosProbability();
  float pionZeroNScore = CVNResults->at(0).GetNpizerosProbability();
  float neutron0Score = CVNResults->at(0).Get0neutronsProbability();
  float neutron1Score = CVNResults->at(0).Get1neutronsProbability();
  float neutron2Score = CVNResults->at(0).Get2neutronsProbability();
  float neutronNScore = CVNResults->at(0).GetNneutronsProbability();

  // Get nu reco information
  const auto numuEOut = e.getValidHandle<dune::EnergyRecoOutput>(fNumuEResultsLabel);
  const auto nueEOut = e.getValidHandle<dune::EnergyRecoOutput>(fNueEResultsLabel);
  const auto NCEOut = e.getValidHandle<dune::EnergyRecoOutput>(fNCEResultsLabel);

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
  fReco.push_back(eventReco);
}

void extrapolation::AddFDData::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  fFile = new HighFive::File(fNDH5FileLoc, HighFive::File::ReadWrite);
}

void extrapolation::AddFDData::endJob()
{
  fFile->createDataSet("fd_reco", fReco);
}

DEFINE_ART_MODULE(extrapolation::AddFDData)
