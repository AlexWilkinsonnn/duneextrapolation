////////////////////////////////////////////////////////////////////////
// Class:       ExportDigits
// Plugin Type: analyzer (Unknown Unknown)
// File:        ExportDigits_module.cc
//
// Generated at Thu Jan 20 08:15:21 2022 by Alexander Wilkinson using cetskelgen
// from  version .
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

#include <memory>
#include <fstream>
#include <iostream>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/raw.h"

namespace extrapolation {
  class ExportDigits;

  struct DepoData;
}

class extrapolation::ExportDigits : public art::EDAnalyzer {
public:
  explicit ExportDigits(fhicl::ParameterSet const& p);

  ExportDigits(ExportDigits const&) = delete;
  ExportDigits(ExportDigits&&) = delete;
  ExportDigits& operator=(ExportDigits const&) = delete;
  ExportDigits& operator=(ExportDigits&&) = delete;

  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  unsigned int fCIndex;
  unsigned int fTIndex;
  unsigned int fPIndex;
  readout::ROPID fRID;
};


extrapolation::ExportDigits::ExportDigits(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fCIndex (p.get<unsigned int>("CryoIndex")),
    fTIndex (p.get<unsigned int>("TpcIndex")),
    fPIndex (p.get<unsigned int>("PlaneIndex"))
{
  consumes<std::vector<raw::RawDigit>>(art::InputTag("tpcrawdecoder", "daq"));
  consumes<std::vector<sim::SimEnergyDeposit>>(art::InputTag("IonAndScint", "EventNumber"));
}

void extrapolation::ExportDigits::analyze(art::Event const& e)
{
  art::Handle<std::vector<raw::RawDigit>> digs;
  e.getByLabel(art::InputTag("tpcrawdecoder", "daq"), digs);
  art::Handle<std::vector<sim::SimEnergyDeposit>> evNumVec;
  e.getByLabel(art::InputTag("IonAndScint", "EventNumber"), evNumVec);

  int evNum = evNumVec->at(0).TrackID();
  std::string fileName = "FD_detsim_" + std::to_string(evNum) + ".csv";

  // Only for collection plane currently
  std::vector<std::vector<short>> ropDigImage(480, std::vector<short>(4492, 0));
  raw::ChannelID_t firstCh = fGeom->FirstChannelInROP(fRID);
  for (const raw::RawDigit& dig : *digs) {
    if (fGeom->ChannelToROP(dig.Channel()) == fRID) {
      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

      // WC detsim still set to 6000 readout window, cannot be changed in fcl (v09_42_00)
      for (unsigned int tick = 0; tick < 4492; tick++) {
        const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;

        ropDigImage[dig.Channel() - firstCh][tick] = adc;
      }
    }
  }

  // Only for collection plane currently
  std::ofstream file;
  file.open(fileName);
  for (int i = 0; i < 480; i++) {
    for (int j = 0; j < 4492; j++) {
      if (ropDigImage[i][j]) {
        file << i << "," << j << "," << ropDigImage[i][j] << "\n";
      }
    }
  }
  file.close();
}

void extrapolation::ExportDigits::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  const geo::CryostatID cID(fCIndex);
  const geo::TPCID tID(cID, fTIndex);
  const geo::PlaneID pID(tID, fPIndex);
  std::cout << pID << "\n";
  fRID = fGeom->WirePlaneToROP(pID);
}

void extrapolation::ExportDigits::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::ExportDigits)
