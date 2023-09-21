////////////////////////////////////////////////////////////////////////
// Class:       AddFDReco
// Plugin Type: analyzer (Unknown Unknown)
// File:        AddFDReco_module.cc
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
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "hep_hpc/hdf5/File.hpp"
#include "hep_hpc/hdf5/make_ntuple.hpp"

#include <string>
#include <vector>
#include <map>

namespace extrapolation {
  class AddFDReco;
}


class extrapolation::AddFDReco : public art::EDAnalyzer {
public:
  explicit AddFDReco(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  // The line "hep_hpc::hdf5::File fFile;" somehow causes the destructor to have noexcept(false)
  // then there is a build error "looser throw specifier ..." coming from the destructor being
  // noexcept in EDAnalyzer. Defining the destructor as below is the only way I found to remove
  // build error. duhduhduhhhhh code????
  virtual ~AddFDReco() noexcept { };

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

  // My function.
  void reset();

private:
  const geo::GeometryCore* fGeom;

  hep_hpc::hdf5::File fFile;

  hep_hpc::hdf5::Ntuple<
    hep_hpc::hdf5::Column<float, 1>, // numu score
    hep_hpc::hdf5::Column<float, 1>, // nue score
    hep_hpc::hdf5::Column<float, 1>, // nc score
    hep_hpc::hdf5::Column<float, 1>, // nutai score
    hep_hpc::hdf5::Column<float, 1>, // antinu score
    hep_hpc::hdf5::Column<float, 1>, // 0 p score
    hep_hpc::hdf5::Column<float, 1>, // 1 p score
    hep_hpc::hdf5::Column<float, 1>, // 2 p score
    hep_hpc::hdf5::Column<float, 1>, // N p score
    hep_hpc::hdf5::Column<float, 1>, // 0 pi score
    hep_hpc::hdf5::Column<float, 1>, // 1 pi score
    hep_hpc::hdf5::Column<float, 1>, // 2 pi score
    hep_hpc::hdf5::Column<float, 1>, // N pi score
    hep_hpc::hdf5::Column<float, 1>, // 0 pi0 score
    hep_hpc::hdf5::Column<float, 1>, // 1 pi0 score
    hep_hpc::hdf5::Column<float, 1>, // 2 pi0 score
    hep_hpc::hdf5::Column<float, 1>, // N pi0 score
    hep_hpc::hdf5::Column<float, 1>, // 0 n score
    hep_hpc::hdf5::Column<float, 1>, // 1 n score
    hep_hpc::hdf5::Column<float, 1>, // 2 n score
    hep_hpc::hdf5::Column<float, 1>, // N n score
    hep_hpc::hdf5::Column<float, 1>, // numu nu E
    hep_hpc::hdf5::Column<float, 1>, // numu had E
    hep_hpc::hdf5::Column<float, 1>, // numu lep E
    hep_hpc::hdf5::Column<int, 1>,   // numu reco method
    hep_hpc::hdf5::Column<int, 1>,   // numu longest track contained
    hep_hpc::hdf5::Column<int, 1>,   // numu longest track momentum method
    hep_hpc::hdf5::Column<float, 1>, // nue nu E
    hep_hpc::hdf5::Column<float, 1>, // nue lep E
    hep_hpc::hdf5::Column<float, 1>, // nue had E
    hep_hpc::hdf5::Column<int, 1>,   // nue reco method
    hep_hpc::hdf5::Column<float, 1>, // nc nu E
    hep_hpc::hdf5::Column<float, 1>, // nc lep E
    hep_hpc::hdf5::Column<float, 1>, // nc had E
    hep_hpc::hdf5::Column<int, 1>    // nc reco method
  >* fEventReco;

  // Labels from fcl
  std::string fCVNResultsLabel;
  std::string fNumuEResultsLabel;
  std::string fNueEResultsLabel;
  std::string fNCEResultsLabel;

  std::string fNDH5FileLoc;
};


extrapolation::AddFDReco::AddFDReco(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fCVNResultsLabel   (p.get<std::string>("CVNResultsLabel")),
    fNumuEResultsLabel (p.get<std::string>("NumuEResultsLabel")),
    fNueEResultsLabel  (p.get<std::string>("NueEResultsLabel")),
    fNCEResultsLabel   (p.get<std::string>("NCEResultsLabel")),
    fNDH5FileLoc       (p.get<std::string>("NDH5FileLoc"))
{
  consumes<std::vector<cvn::Result>>(fCVNResultsLabel);

  consumes<dune::EnergyRecoOutput>(fNumuEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNueEResultsLabel);
  consumes<dune::EnergyRecoOutput>(fNCEResultsLabel);

  consumes<std::vector<recob::Hit>>(fNCEResultsLabel);
}

void extrapolation::AddFDReco::analyze(art::Event const& e)
{
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

  fEventReco->insert(
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
  );
}

void extrapolation::AddFDReco::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  fFile = hep_hpc::hdf5::File(fNDH5FileLoc, H5F_ACC_RDWR);

  fEventReco = new hep_hpc::hdf5::Ntuple(
    hep_hpc::hdf5::make_ntuple(
      {fFile, "fd_reco", 1000},
      hep_hpc::hdf5::make_scalar_column<float>("numu_score"),
      hep_hpc::hdf5::make_scalar_column<float>("nue_score"),
      hep_hpc::hdf5::make_scalar_column<float>("nc_score"),
      hep_hpc::hdf5::make_scalar_column<float>("nutau_score"),
      hep_hpc::hdf5::make_scalar_column<float>("antinu_score"),
      hep_hpc::hdf5::make_scalar_column<float>("0_p_score"),
      hep_hpc::hdf5::make_scalar_column<float>("1_p_score"),
      hep_hpc::hdf5::make_scalar_column<float>("2_p_score"),
      hep_hpc::hdf5::make_scalar_column<float>("N_p_score"),
      hep_hpc::hdf5::make_scalar_column<float>("0_pi_score"),
      hep_hpc::hdf5::make_scalar_column<float>("1_pi_score"),
      hep_hpc::hdf5::make_scalar_column<float>("2_pi_score"),
      hep_hpc::hdf5::make_scalar_column<float>("N_pi_score"),
      hep_hpc::hdf5::make_scalar_column<float>("0_pi0_score"),
      hep_hpc::hdf5::make_scalar_column<float>("1_pi0_score"),
      hep_hpc::hdf5::make_scalar_column<float>("2_pi0_score"),
      hep_hpc::hdf5::make_scalar_column<float>("N_pi0_score"),
      hep_hpc::hdf5::make_scalar_column<float>("0_n_score"),
      hep_hpc::hdf5::make_scalar_column<float>("1_n_score"),
      hep_hpc::hdf5::make_scalar_column<float>("2_n_score"),
      hep_hpc::hdf5::make_scalar_column<float>("N_n_score"),
      hep_hpc::hdf5::make_scalar_column<float>("numu_nu_E"),
      hep_hpc::hdf5::make_scalar_column<float>("numu_had_E"),
      hep_hpc::hdf5::make_scalar_column<float>("numu_lep_E"),
      hep_hpc::hdf5::make_scalar_column<int>("numu_reco_method"),
      hep_hpc::hdf5::make_scalar_column<int>("numu_longest_track_contained"),
      hep_hpc::hdf5::make_scalar_column<int>("numu_longest_track_mom_method"),
      hep_hpc::hdf5::make_scalar_column<float>("nue_nu_E"),
      hep_hpc::hdf5::make_scalar_column<float>("nue_had_E"),
      hep_hpc::hdf5::make_scalar_column<float>("nue_lep_E"),
      hep_hpc::hdf5::make_scalar_column<int>("nue_reco_method"),
      hep_hpc::hdf5::make_scalar_column<float>("nc_nu_E"),
      hep_hpc::hdf5::make_scalar_column<float>("nc_had_E"),
      hep_hpc::hdf5::make_scalar_column<float>("nc_lep_E"),
      hep_hpc::hdf5::make_scalar_column<int>("nc_reco_method")
    )
  );

}

void extrapolation::AddFDReco::endJob()
{
  delete fEventReco;
  fFile.close();
}

DEFINE_ART_MODULE(extrapolation::AddFDReco)
