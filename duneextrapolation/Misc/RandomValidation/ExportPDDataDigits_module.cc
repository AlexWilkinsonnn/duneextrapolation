////////////////////////////////////////////////////////////////////////
// Class:       ExportPDDataDigits
// Plugin Type: analyzer (Unknown Unknown)
// File:        ExportPDDataDigits_module.cc
//
// Generated Wed Aug 21 2024
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
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace extrapolation {
  class ExportPDDataDigits;
}

class extrapolation::ExportPDDataDigits : public art::EDAnalyzer {
public:
  explicit ExportPDDataDigits(fhicl::ParameterSet const& p);

  ExportPDDataDigits(ExportPDDataDigits const&) = delete;
  ExportPDDataDigits(ExportPDDataDigits&&) = delete;
  ExportPDDataDigits& operator=(ExportPDDataDigits const&) = delete;
  ExportPDDataDigits& operator=(ExportPDDataDigits&&) = delete;

  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

  void reset();
private:
  const geo::GeometryCore* fGeom;

  TTree*                        fTreeDigits;
  // Need to use int instead of short to avoid providing a definition file for ROOT
  std::vector<std::vector<int>> fDigitsZ;
  std::vector<std::vector<int>> fDigitsU;
  std::vector<std::vector<int>> fDigitsV;
  std::vector<int>              fPedsZ;
  std::vector<int>              fPedsU;
  std::vector<int>              fPedsV;
  int                           fEventNum;

  unsigned int fCIndex;
  unsigned int fTIndex;
  geo::PlaneID   fPIDZ;
  geo::PlaneID   fPIDU;
  geo::PlaneID   fPIDV;
  readout::ROPID fRIDZ;
  readout::ROPID fRIDU;
  readout::ROPID fRIDV;
};

extrapolation::ExportPDDataDigits::ExportPDDataDigits(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fCIndex     (p.get<unsigned int>("CryoIndex")),
    fTIndex     (p.get<unsigned int>("TpcIndex"))
{
  consumes<std::vector<raw::RawDigit>>(art::InputTag("tpcrawdecoder", "daq"));

  art::ServiceHandle<art::TFileService> tfs;

  fTreeDigits = tfs->make<TTree>("digits", "digits");
  fTreeDigits->Branch("digit_vecsZ", &fDigitsZ);
  fTreeDigits->Branch("digit_vecsU", &fDigitsU);
  fTreeDigits->Branch("digit_vecsV", &fDigitsV);
  fTreeDigits->Branch("digit_pedZ", &fPedsZ);
  fTreeDigits->Branch("digit_pedU", &fPedsU);
  fTreeDigits->Branch("digit_pedV", &fPedsV);
  fTreeDigits->Branch("eventid", &fEventNum);
}

void extrapolation::ExportPDDataDigits::analyze(art::Event const& e)
{
  this->reset();

  art::Handle<std::vector<raw::RawDigit>> digs;
  e.getByLabel(art::InputTag("tpcrawdecoder", "daq"), digs);

  std::cout << e.run() << " - " << e.subRun() << " - " << e.event() << "\n";
  fEventNum = (int)e.event();

  fDigitsZ = std::vector<std::vector<int>>(480);
  fDigitsU = std::vector<std::vector<int>>(800);
  fDigitsV = std::vector<std::vector<int>>(800);
  fPedsZ = std::vector<int>(480);
  fPedsU = std::vector<int>(800);
  fPedsV = std::vector<int>(800);
  raw::ChannelID_t firstChZ = fGeom->FirstChannelInROP(fRIDZ);
  raw::ChannelID_t firstChU = fGeom->FirstChannelInROP(fRIDU);
  raw::ChannelID_t firstChV = fGeom->FirstChannelInROP(fRIDV);
  for (const raw::RawDigit& dig : *digs) {
    if (fGeom->ChannelToROP(dig.Channel()) == fRIDZ) {
      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

      // WC detsim still set to 6000 readout window, cannot be changed in fcl (v09_42_00)
      fDigitsZ[dig.Channel() - firstChZ] = std::vector<int>(4492);
      fPedsZ[dig.Channel() - firstChZ] = (int)dig.GetPedestal();
      for (unsigned int tick = 0; tick < 4492; tick++) {
        // const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;
        const short adc = adcs[tick];

        fDigitsZ[dig.Channel() - firstChZ][tick] = (int)adc;
      }
    }
    else if (fGeom->ChannelToROP(dig.Channel()) == fRIDU) {
      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

      // WC detsim still set to 6000 readout window, cannot be changed in fcl (v09_42_00)
      fDigitsU[dig.Channel() - firstChU] = std::vector<int>(4492);
      fPedsU[dig.Channel() - firstChU] = (int)dig.GetPedestal();
      for (unsigned int tick = 0; tick < 4492; tick++) {
        // const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;
        const short adc = adcs[tick];

        fDigitsU[dig.Channel() - firstChU][tick] = (int)adc;
      }
    }
    else if (fGeom->ChannelToROP(dig.Channel()) == fRIDV) {
      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

      // WC detsim still set to 6000 readout window, cannot be changed in fcl (v09_42_00)
      fDigitsV[dig.Channel() - firstChV] = std::vector<int>(4492);
      fPedsV[dig.Channel() - firstChV] = (int)dig.GetPedestal();
      for (unsigned int tick = 0; tick < 4492; tick++) {
        // const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;
        const short adc = adcs[tick];

        fDigitsV[dig.Channel() - firstChV][tick] = (int)adc;
      }
    }
  }

  fTreeDigits->Fill();
}

void extrapolation::ExportPDDataDigits::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  const geo::CryostatID cID(fCIndex);
  const geo::TPCID tID(cID, fTIndex);
  for (geo::PlaneID pID : fGeom->Iterate<geo::PlaneID>(tID)) {
    if (fGeom->View(pID) == geo::kZ) {
      std::cout << "Z plane: " << pID << "\n";
      fPIDZ = pID;
      fRIDZ = fGeom->WirePlaneToROP(pID);
    }
    else if (fGeom->View(pID) == geo::kU) {
      std::cout << "U plane: " << pID << "\n";
      fPIDU = pID;
      fRIDU = fGeom->WirePlaneToROP(pID);
    }
    else if (fGeom->View(pID) == geo::kV) {
      std::cout << "V plane: " << pID << "\n";
      fPIDV = pID;
      fRIDV = fGeom->WirePlaneToROP(pID);
    }
  }
}

void extrapolation::ExportPDDataDigits::endJob()
{
}

void extrapolation::ExportPDDataDigits::reset()
{
  fDigitsZ.clear();
  fDigitsU.clear();
  fDigitsV.clear();
  fPedsZ.clear();
  fPedsU.clear();
  fPedsV.clear();
  fEventNum = -1;
}

DEFINE_ART_MODULE(extrapolation::ExportPDDataDigits)
