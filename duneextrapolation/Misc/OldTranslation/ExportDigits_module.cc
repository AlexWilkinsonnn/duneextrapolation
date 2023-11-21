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
  class ExportDigits;
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

  void reset();
private:
  const geo::GeometryCore* fGeom;

  TTree*                        fTreeDigits;
  // Need to use int instead of short to avoid providing a definition file for ROOT
  std::vector<std::vector<int>> fDigitsZ;
  std::vector<std::vector<int>> fDigitsU;
  std::vector<std::vector<int>> fDigitsV;
  int                           fEventNum;

  TTree*                        fTreeSEDs;
  std::vector<std::vector<int>> fSEDsZ;
  std::vector<std::vector<int>> fSEDsU;
  std::vector<std::vector<int>> fSEDsV;


  unsigned int fCIndex;
  unsigned int fTIndex;
  geo::PlaneID   fPIDZ;
  geo::PlaneID   fPIDU;
  geo::PlaneID   fPIDV;
  readout::ROPID fRIDZ;
  readout::ROPID fRIDU;
  readout::ROPID fRIDV;

  bool fExportSEDs;
};


extrapolation::ExportDigits::ExportDigits(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fCIndex     (p.get<unsigned int>("CryoIndex")),
    fTIndex     (p.get<unsigned int>("TpcIndex")),
    fExportSEDs (p.get<bool>("ExportSEDs"))
{
  consumes<std::vector<raw::RawDigit>>(art::InputTag("tpcrawdecoder", "daq"));
  consumes<std::vector<sim::SimEnergyDeposit>>(art::InputTag("IonAndScint", "EventNumber"));

  art::ServiceHandle<art::TFileService> tfs;

  fTreeDigits = tfs->make<TTree>("digits", "digits");
  fTreeDigits->Branch("digit_vecsZ", &fDigitsZ);
  fTreeDigits->Branch("digit_vecsU", &fDigitsU);
  fTreeDigits->Branch("digit_vecsV", &fDigitsV);
  fTreeDigits->Branch("eventid", &fEventNum);

  if (fExportSEDs) {
    consumes<std::vector<sim::SimEnergyDeposit>>(art::InputTag("IonAndScint", ""));

    fTreeSEDs = tfs->make<TTree>("seds", "seds");
    fTreeSEDs->Branch("sedsZ", &fSEDsZ);
    fTreeSEDs->Branch("sedsU", &fSEDsU);
    fTreeSEDs->Branch("sedsV", &fSEDsV);
    fTreeSEDs->Branch("eventid", &fEventNum);
  }
}

void extrapolation::ExportDigits::analyze(art::Event const& e)
{
  this->reset();

  art::Handle<std::vector<raw::RawDigit>> digs;
  e.getByLabel(art::InputTag("tpcrawdecoder", "daq"), digs);
  art::Handle<std::vector<sim::SimEnergyDeposit>> evNumVec;
  e.getByLabel(art::InputTag("IonAndScint", "EventNumber"), evNumVec);

  int evNum = evNumVec->at(0).TrackID();
  fEventNum = evNum;

  fDigitsZ = std::vector<std::vector<int>>(480);
  fDigitsU = std::vector<std::vector<int>>(800);
  fDigitsV = std::vector<std::vector<int>>(800);
  raw::ChannelID_t firstChZ = fGeom->FirstChannelInROP(fRIDZ);
  raw::ChannelID_t firstChU = fGeom->FirstChannelInROP(fRIDU);
  raw::ChannelID_t firstChV = fGeom->FirstChannelInROP(fRIDV);
  for (const raw::RawDigit& dig : *digs) {
    if (fGeom->ChannelToROP(dig.Channel()) == fRIDZ) {
      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

      // WC detsim still set to 6000 readout window, cannot be changed in fcl (v09_42_00)
      fDigitsZ[dig.Channel() - firstChZ] = std::vector<int>(4492);
      for (unsigned int tick = 0; tick < 4492; tick++) {
        const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;

        fDigitsZ[dig.Channel() - firstChZ][tick] = (int)adc;
      }
    }
    else if (fGeom->ChannelToROP(dig.Channel()) == fRIDU) {
      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

      // WC detsim still set to 6000 readout window, cannot be changed in fcl (v09_42_00)
      fDigitsU[dig.Channel() - firstChU] = std::vector<int>(4492);
      for (unsigned int tick = 0; tick < 4492; tick++) {
        const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;

        fDigitsU[dig.Channel() - firstChU][tick] = (int)adc;
      }
    }
    else if (fGeom->ChannelToROP(dig.Channel()) == fRIDV) {
      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

      // WC detsim still set to 6000 readout window, cannot be changed in fcl (v09_42_00)
      fDigitsV[dig.Channel() - firstChV] = std::vector<int>(4492);
      for (unsigned int tick = 0; tick < 4492; tick++) {
        const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;

        fDigitsV[dig.Channel() - firstChV][tick] = (int)adc;
      }
    }
  }

  fTreeDigits->Fill();

  if (fExportSEDs) {
    const auto SEDs = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(art::InputTag("IonAndScint", ""));

    for (auto& SED : *SEDs) {
      auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

      raw::ChannelID_t chZ = fGeom->NearestChannel(SED.MidPoint(), fPIDZ);
      raw::ChannelID_t chU = fGeom->NearestChannel(SED.MidPoint(), fPIDU);
      raw::ChannelID_t chV = fGeom->NearestChannel(SED.MidPoint(), fPIDV);

      int chLocalZ = (int)(chZ - firstChZ);
      int chLocalU = (int)(chU - firstChU);
      int chLocalV = (int)(chV - firstChV);

      double tickRawZ = detProp.ConvertXToTicks(SED.MidPoint().X(), fPIDZ) + clockData.TPCG4Time2TDC(SED.Time());
      double tickRawU = detProp.ConvertXToTicks(SED.MidPoint().X(), fPIDU) + clockData.TPCG4Time2TDC(SED.Time());
      double tickRawV = detProp.ConvertXToTicks(SED.MidPoint().X(), fPIDV) + clockData.TPCG4Time2TDC(SED.Time());

      tickRawZ -= 7.8;
      tickRawU -= 10.1;
      tickRawV -= 10.9;

      int tickZ = (int)tickRawZ;
      int tickU = (int)tickRawU;
      int tickV = (int)tickRawV;

      fSEDsZ.push_back({chLocalZ, tickZ, SED.NumElectrons()});
      fSEDsU.push_back({chLocalU, tickU, SED.NumElectrons()});
      fSEDsV.push_back({chLocalV, tickV, SED.NumElectrons()});
    }

    fTreeSEDs->Fill();
  }
}

void extrapolation::ExportDigits::beginJob()
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

void extrapolation::ExportDigits::endJob()
{
}

void extrapolation::ExportDigits::reset()
{
  fDigitsZ.clear();
  fDigitsU.clear();
  fDigitsV.clear();
  fEventNum = -1;

  if (fExportSEDs) {
    fSEDsZ.clear();
    fSEDsU.clear();
    fSEDsV.clear();
  }
}

DEFINE_ART_MODULE(extrapolation::ExportDigits)

