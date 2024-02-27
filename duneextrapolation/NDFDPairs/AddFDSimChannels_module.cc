////////////////////////////////////////////////////////////////////////
// Class:       AddFDSimChannels
// Plugin Type: analyzer (Unknown Unknown)
// File:        AddFDSimChannels.cc
//
// Created on 27 Feb 24 Alex Wilkinson
// Dump SimChannels to the ND-FD pair HDF5 file.
// Used to pass SimChannels back to old dunetpc for detsim+reco to get
// TDR era reco-reco pairs. Cannot do SimEnergyDeposit->SimChannel in
// the old dunetpc since it is before the LArG4 refactor
// (not impossible just far too much work).
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

#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "highfive/H5DataSet.hpp"
#include "highfive/H5File.hpp"
#include "highfive/H5Group.hpp"
#include "highfive/H5Object.hpp"
#include "highfive/H5DataType.hpp"
#include "highfive/H5DataSpace.hpp"

#include <string>
#include <vector>
#include <map>

typedef struct simChannelIonisationFD {
  int trackID;
  unsigned int tdc;
  double numberElectrons;
  double x;
  double y;
  double z;
  double energy;
} simChannelIonisationFD;

HighFive::CompoundType make_simChannelIonisationFD() {
  return {
    {"trackID", HighFive::AtomicType<int>{}},
    {"tdc", HighFive::AtomicType<unsigned int>{}},
    {"numberElectrons", HighFive::AtomicType<double>{}},
    {"x", HighFive::AtomicType<double>{}},
    {"y", HighFive::AtomicType<double>{}},
    {"z", HighFive::AtomicType<double>{}},
    {"energy", HighFive::AtomicType<double>{}}
  };
}

HIGHFIVE_REGISTER_TYPE(simChannelIonisationFD, make_simChannelIonisationFD)

namespace extrapolation {
  class AddFDSimChannels;
}

class extrapolation::AddFDSimChannels : public art::EDAnalyzer {
public:
  explicit AddFDSimChannels(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  AddFDSimChannels(AddFDSimChannels const&) = delete;
  AddFDSimChannels(AddFDSimChannels&&) = delete;
  AddFDSimChannels& operator=(AddFDSimChannels const&) = delete;
  AddFDSimChannels& operator=(AddFDSimChannels&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  // Methods
  
  // Members
  HighFive::File* fFile;

  // Product labels
  std::string fEventIDSEDLabel;
  std::string fSimChannelsLabel;

  // fcl params
  std::string fNDFDH5FileLoc;
};


extrapolation::AddFDSimChannels::AddFDSimChannels(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fEventIDSEDLabel  (p.get<std::string>("EventIDSEDLabel")),
    fSimChannelsLabel (p.get<std::string>("CVNResultsLabel")),
    fNDFDH5FileLoc    (p.get<std::string>("NDFDH5FileLoc"))
{
  consumes<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);
  consumes<std::vector<sim::SimChannel>>(fSimChannelsLabel);
}

void extrapolation::AddFDSimChannels::analyze(art::Event const& e)
{
  // Get eventID
  const auto eventIDSED = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);
  int eventID = (*eventIDSED)[0].TrackID();

  const auto SCs = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelsLabel);

  for (const sim::SimChannel SC : *SCs) {
    const raw::ChannelID_t ch = SC.Channel();

    std::vector<simChannelIonisationFD> ionisations;
    for (const sim::TDCIDE ionisationTDCIDEs : SC.TDCIDEMap()) {
      unsigned int tdc = (unsigned int)ionisationTDCIDEs.first;

      for (const sim::IDE ionisationIDE : ionisationTDCIDEs.second) {
        const simChannelIonisationFD ionisation = {
          ionisationIDE.trackID,
          tdc,
          (double)ionisationIDE.numElectrons,
          (double)ionisationIDE.x,
          (double)ionisationIDE.y,
          (double)ionisationIDE.z,
          (double)ionisationIDE.energy
        };
        ionisations.push_back(ionisation);
      }
    }

    const std::string groupPath = 
      "fd_simchannels/" + std::to_string(eventID) + "/" + std::to_string(ch);

    fFile->createDataSet(groupPath, ionisations);
  }
}

void extrapolation::AddFDSimChannels::beginJob()
{
  fFile = new HighFive::File(fNDFDH5FileLoc, HighFive::File::ReadWrite);
}

DEFINE_ART_MODULE(extrapolation::AddFDSimChannels)
