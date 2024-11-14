////////////////////////////////////////////////////////////////////////
// Class:       GetNoise
// Plugin Type: analyzer (Unknown Unknown)
// File:        GetNoise_module.cc
//
// Crated on 14 Nov 24 Alex Wilkinson
// Gets one plane of only noise for U,V,Z and writes to a h5 file.
// This is used to add fake noise to the predicted FD response loaded
// in with LoadFDResp.
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
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/raw.h"

#include "highfive/H5DataSet.hpp"
#include "highfive/H5File.hpp"
#include "highfive/H5Group.hpp"
#include "highfive/H5Object.hpp"
#include "highfive/H5DataType.hpp"
#include "highfive/H5DataSpace.hpp"

#include <string>
#include <vector>
#include <set>
#include <iomanip>

namespace extrapolation {
  class GetNoise;
}

class extrapolation::GetNoise : public art::EDAnalyzer {
public:
  explicit GetNoise(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  GetNoise(GetNoise const&) = delete;
  GetNoise(GetNoise&&) = delete;
  GetNoise& operator=(GetNoise const&) = delete;
  GetNoise& operator=(GetNoise&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Methods
  bool inWireCellBoundingBox(const double x, const double y, const double z);
  double roundDoubleSigFigs(const double val, const int nSigfigs);

  // Members
  const geo::GeometryCore* fGeom;

  HighFive::File* fFile;

  // Product labels
  std::string fFDSEDLabel;
  std::string fRawDigitLabel;

  // See config fcl for description of these members
  std::string fNDFDH5FileName;
  std::vector<std::vector<std::vector<double>>> fWireCellAPABoundingBoxes;
  unsigned int fMaxTick;
};


extrapolation::GetNoise::GetNoise(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fFDSEDLabel                (p.get<std::string>("FDSEDLabel")),
    fRawDigitLabel             (p.get<std::string>("RawDigitLabel")),
    fNDFDH5FileName            (p.get<std::string>("NDFDH5FileName")),
    fWireCellAPABoundingBoxes  (p.get<std::vector<std::vector<std::vector<double>>>>("WireCellAPABoundingBoxes")),
    fMaxTick                   (p.get<unsigned int>("MaxTick"))
{
  consumes<std::vector<raw::RawDigit>>(fRawDigitLabel);
  consumes<std::vector<sim::SimEnergyDeposit>>(fFDSEDLabel);

  fMaxTick = std::min(fMaxTick, (unsigned int)6000);
}

void extrapolation::GetNoise::analyze(art::Event const& e)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
  
  // Find ROPIDs that have signal
  std::set<readout::ROPID> activeRIDs;
  const auto SEDs = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fFDSEDLabel);
  for (const auto& SED : *SEDs) { 
    if (!inWireCellBoundingBox(SED.X(), SED.Y(), SED.Z())) {
      continue;
    }
    const geo::Point_t SEDLoc = SED.MidPoint();
    const geo::TPCID tID = fGeom->PositionToTPCID(SEDLoc);
    for (const geo::PlaneID pID : fGeom->Iterate<geo::PlaneID>(tID)) {
      const readout::ROPID rID = fGeom->WirePlaneToROP(pID);
      activeRIDs.insert(rID);
    }
  }

  // Get the noise-only RawDigits for each plane
  const auto digits = e.getValidHandle<std::vector<raw::RawDigit>>(fRawDigitLabel);
  std::map<readout::ROPID, std::vector<std::vector<short>>> rIDDigs;
  for (const raw::RawDigit& dig : *digits) { 
    const readout::ROPID rID = fGeom->ChannelToROP(dig.Channel());
    if (activeRIDs.find(rID) != activeRIDs.end()) {
      continue;
    }

    if (rIDDigs.find(rID) == rIDDigs.end()) {
      rIDDigs[rID] = 
        std::vector<std::vector<short>>(fGeom->Nchannels(rID), std::vector<short>(fMaxTick, 0));
    }

    raw::RawDigit::ADCvector_t adcs(dig.Samples());
    raw::Uncompress(dig.ADCs(), adcs, dig.Compression());
    for (unsigned int tick = 0; tick < fMaxTick; tick++) { 
      const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;
      rIDDigs[rID][dig.Channel() - fGeom->FirstChannelInROP(rID)][tick] = adc;
    }
  }

  // Write out one example of noise-only for each plane
  std::set<geo::View_t> doneViews;
  for (const auto& rID_adcs : rIDDigs) {
    const readout::ROPID rID = rID_adcs.first;
    if (doneViews.find(fGeom->View(rID)) != doneViews.end()){
      continue;
    }

    const std::vector<std::vector<short>> adcs = rID_adcs.second;
    std::string groupPath = "noises/";
    if (fGeom->View(rID) == geo::kZ) {
      groupPath = groupPath + "Z";
    }
    else if (fGeom->View(rID) == geo::kU) {
      groupPath = groupPath + "U";
    }
    else if (fGeom->View(rID) == geo::kV) {
      groupPath = groupPath + "V";
    }
    fFile->createDataSet(groupPath, adcs);
    doneViews.insert(fGeom->View(rID));
  }
}

void extrapolation::GetNoise::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  fFile = new HighFive::File(fNDFDH5FileName, HighFive::File::Create);
}

void extrapolation::GetNoise::endJob()
{
}

bool extrapolation::GetNoise::inWireCellBoundingBox(
  const double x, const double y, const double z
)
{
  // WireCell bounding box is only known to 6 significant figures
  const double xRounded = roundDoubleSigFigs(x, 6);
  const double yRounded = roundDoubleSigFigs(y, 6);
  const double zRounded = roundDoubleSigFigs(z, 6);
  for (const std::vector<std::vector<double>>& range : fWireCellAPABoundingBoxes) {
    if (
      xRounded > range[0][0] && xRounded < range[1][0] &&
      yRounded > range[0][1] && yRounded < range[1][1] &&
      zRounded > range[0][2] && zRounded < range[1][2]
    ) {
      return true;
    }
  }
  return false;
}

double extrapolation::GetNoise::roundDoubleSigFigs(const double val, const int nSigFigs)
{
  // *yoink* (https://cplusplus.com/forum/beginner/274016/)
  std::stringstream ss;
  ss << std::scientific << std::setprecision(nSigFigs - 1) << val;
  return stod(ss.str());
}

DEFINE_ART_MODULE(extrapolation::GetNoise)
