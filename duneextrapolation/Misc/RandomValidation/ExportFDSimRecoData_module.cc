////////////////////////////////////////////////////////////////////////
// Class:       ExportFDSimRecoData
// Plugin Type: analyzer (Unknown Unknown)
// File:        ExportFDSimRecoData_module.cc
//
// Generated Wed Aug 23 2024
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

// art
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

// LArSoft
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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// STL
#include <memory>
#include <fstream>
#include <iostream>
#include <vector>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace extrapolation {
  class ExportFDSimRecoData;
}

class extrapolation::ExportFDSimRecoData : public art::EDAnalyzer {
public:
  explicit ExportFDSimRecoData(fhicl::ParameterSet const& p);

  ExportFDSimRecoData(ExportFDSimRecoData const&) = delete;
  ExportFDSimRecoData(ExportFDSimRecoData&&) = delete;
  ExportFDSimRecoData& operator=(ExportFDSimRecoData const&) = delete;
  ExportFDSimRecoData& operator=(ExportFDSimRecoData&&) = delete;

  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

private:
  void reset();
  void getVtxROPs(
    const art::Event &e, readout::ROPID &rIDZ, readout::ROPID &rIDU, readout::ROPID &rIDV
  );
  void fillDigits(const raw::RawDigit &dig, std::vector<int> &digVec);

  const geo::GeometryCore* fGeom;

  bool fExportDigits;

  TTree*                        fTreeSimReco;
  // Need to use int instead of short to avoid providing a definition file for ROOT
  std::vector<std::vector<int>> fDigitsZ;
  std::vector<std::vector<int>> fDigitsU;
  std::vector<std::vector<int>> fDigitsV;
  int                           fEventNum;
};

extrapolation::ExportFDSimRecoData::ExportFDSimRecoData(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fExportDigits (p.get<bool>("ExportDigits"))
{
  consumes<std::vector<simb::MCTruth>>(art::InputTag("generator"));

  art::ServiceHandle<art::TFileService> tfs;

  fTreeSimReco = tfs->make<TTree>("simreco", "simreco");

  if (fExportDigits) {
    consumes<std::vector<raw::RawDigit>>(art::InputTag("tpcrawdecoder", "daq"));
    fTreeSimReco->Branch("digit_vecsZ", &fDigitsZ);
    fTreeSimReco->Branch("digit_vecsU", &fDigitsU);
    fTreeSimReco->Branch("digit_vecsV", &fDigitsV);
    fTreeSimReco->Branch("eventid", &fEventNum);
  }
}

void extrapolation::ExportFDSimRecoData::analyze(art::Event const& e)
{
  this->reset();

  std::cout << e.run() << " - " << e.subRun() << " - " << e.event() << "\n";
  fEventNum = (int)e.event();

  readout::ROPID rIDZ, rIDU, rIDV;
  getVtxROPs(e, rIDZ, rIDU, rIDV);

  if (!fGeom->HasROP(rIDZ) || !fGeom->HasROP(rIDU) || !fGeom->HasROP(rIDV)) {
    return;
  }
  
  if (fExportDigits) {
    art::Handle<std::vector<raw::RawDigit>> digs;
    e.getByLabel(art::InputTag("tpcrawdecoder", "daq"), digs);

    fDigitsZ = std::vector<std::vector<int>>(480);
    fDigitsU = std::vector<std::vector<int>>(800);
    fDigitsV = std::vector<std::vector<int>>(800);
    raw::ChannelID_t firstChZ = fGeom->FirstChannelInROP(rIDZ);
    raw::ChannelID_t firstChU = fGeom->FirstChannelInROP(rIDU);
    raw::ChannelID_t firstChV = fGeom->FirstChannelInROP(rIDV);
    for (const raw::RawDigit& dig : *digs) {
      if (fGeom->ChannelToROP(dig.Channel()) == rIDZ) {
        // WC detsim still set to 6000 readout window, cannot be changed in fcl (v09_42_00)
        fDigitsZ[dig.Channel() - firstChZ] = std::vector<int>(4492);
        fillDigits(dig, fDigitsZ[dig.Channel() - firstChZ]);
      }
      else if (fGeom->ChannelToROP(dig.Channel()) == rIDU) {
        fDigitsU[dig.Channel() - firstChU] = std::vector<int>(4492);
        fillDigits(dig, fDigitsU[dig.Channel() - firstChU]);
      }
      else if (fGeom->ChannelToROP(dig.Channel()) == rIDV) {
        fDigitsV[dig.Channel() - firstChV] = std::vector<int>(4492);
        fillDigits(dig, fDigitsV[dig.Channel() - firstChV]);
      }
    }
  }

  fTreeSimReco->Fill();
}

void extrapolation::ExportFDSimRecoData::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
}

void extrapolation::ExportFDSimRecoData::endJob()
{
}

void extrapolation::ExportFDSimRecoData::reset()
{
  fDigitsZ.clear();
  fDigitsU.clear();
  fDigitsV.clear();
  fEventNum = -1;
}

void extrapolation::ExportFDSimRecoData::getVtxROPs(
  const art::Event &e, readout::ROPID &rIDZ, readout::ROPID &rIDU, readout::ROPID &rIDV
)
{
  art::Handle<std::vector<simb::MCTruth>> mcTruth;
  e.getByLabel(art::InputTag("generator"), mcTruth);
  const simb::MCParticle &mcNuPart = mcTruth->at(0).GetNeutrino().Nu();
  const geo::Point_t vtx(mcNuPart.Vx(), mcNuPart.Vy(), mcNuPart.Vz());
  const geo::TPCID vtxTID = fGeom->PositionToTPCID(vtx);
  std::cout << "Nu vertex is (" << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() << ")\n";

  for (const geo::PlaneID &pID : fGeom->Iterate<geo::PlaneID>(vtxTID)) {
    if (fGeom->View(pID) == geo::kZ) {
      std::cout << "Z plane: " << pID << "\n";
      rIDZ = fGeom->WirePlaneToROP(pID);
    }
    else if (fGeom->View(pID) == geo::kU) {
      std::cout << "U plane: " << pID << "\n";
      rIDU = fGeom->WirePlaneToROP(pID);
    }
    else if (fGeom->View(pID) == geo::kV) {
      std::cout << "V plane: " << pID << "\n";
      rIDV = fGeom->WirePlaneToROP(pID);
    }
  }
}

void extrapolation::ExportFDSimRecoData::fillDigits(
  const raw::RawDigit &dig, std::vector<int> &digVec
)
{
  raw::RawDigit::ADCvector_t adcs(dig.Samples());
  raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

  for (unsigned int tick = 0; tick < digVec.size(); tick++) {
    const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;
    digVec[tick] = (int)adc;
  }
}

DEFINE_ART_MODULE(extrapolation::ExportFDSimRecoData)
