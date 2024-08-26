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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// STL
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
  void fillWires(const recob::Wire &wire, std::vector<float> &wireVec);

  const geo::GeometryCore* fGeom;

  std::string fSPLabel;

  bool fExportDigits;
  bool fExportSPWires;
  bool fExportHits;

  TTree*                        fTreeSimReco;
  // Need to use int instead of short to avoid providing a definition file for ROOT
  std::vector<std::vector<int>>   fDigitsZ;
  std::vector<std::vector<int>>   fDigitsU;
  std::vector<std::vector<int>>   fDigitsV;
  std::vector<std::vector<float>> fWiresZ;
  std::vector<std::vector<float>> fWiresU;
  std::vector<std::vector<float>> fWiresV;
  std::vector<int>                fHitROPs;
  std::vector<int>                fHitChs;
  std::vector<float>              fHitPeakTimes;
  std::vector<int>                fHitStartTicks;
  std::vector<int>                fHitEndTicks;
  std::vector<float>              fHitPeakSigmas;
  std::vector<float>              fHitRMSs;
  std::vector<float>              fHitPeakAmplitudes;
  int                             fEventNum;
};

extrapolation::ExportFDSimRecoData::ExportFDSimRecoData(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSPLabel       (p.get<std::string>("SPLabel")),
    fExportDigits  (p.get<bool>("ExportDigits")),
    fExportSPWires (p.get<bool>("ExportSPWires")),
    fExportHits    (p.get<bool>("ExportHits"))
{
  consumes<std::vector<simb::MCTruth>>(art::InputTag("generator"));

  art::ServiceHandle<art::TFileService> tfs;

  fTreeSimReco = tfs->make<TTree>("simreco", "simreco");
  fTreeSimReco->Branch("eventid", &fEventNum);

  if (fExportDigits) {
    consumes<std::vector<raw::RawDigit>>(art::InputTag("tpcrawdecoder", "daq"));
    fTreeSimReco->Branch("digit_vecsZ", &fDigitsZ);
    fTreeSimReco->Branch("digit_vecsU", &fDigitsU);
    fTreeSimReco->Branch("digit_vecsV", &fDigitsV);
  }

  if (fExportSPWires) {
    consumes<std::vector<recob::Wire>>(fSPLabel);
    fTreeSimReco->Branch("wire_vecsZ", &fWiresZ);
    fTreeSimReco->Branch("wire_vecsU", &fWiresU);
    fTreeSimReco->Branch("wire_vecsV", &fWiresV);
  }

  if (fExportHits) {
    consumes<std::vector<recob::Hit>>(art::InputTag("gaushit"));
    fTreeSimReco->Branch("hit_planes", &fHitROPs);
    fTreeSimReco->Branch("hit_chs", &fHitChs);
    fTreeSimReco->Branch("hit_peaktimes", &fHitPeakTimes);
    fTreeSimReco->Branch("hit_peaktimes", &fHitPeakTimes);
    fTreeSimReco->Branch("hit_startticks", &fHitStartTicks);
    fTreeSimReco->Branch("hit_endticks", &fHitEndTicks);
    fTreeSimReco->Branch("hit_peaksigmas", &fHitPeakSigmas);
    fTreeSimReco->Branch("hit_rmss", &fHitRMSs);
    fTreeSimReco->Branch("hit_peakamplitudes", &fHitPeakAmplitudes);
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
  raw::ChannelID_t firstChZ = fGeom->FirstChannelInROP(rIDZ);
  raw::ChannelID_t firstChU = fGeom->FirstChannelInROP(rIDU);
  raw::ChannelID_t firstChV = fGeom->FirstChannelInROP(rIDV);
  
  if (fExportDigits) {
    art::Handle<std::vector<raw::RawDigit>> digs;
    e.getByLabel(art::InputTag("tpcrawdecoder", "daq"), digs);

    fDigitsZ = std::vector<std::vector<int>>(480);
    fDigitsU = std::vector<std::vector<int>>(800);
    fDigitsV = std::vector<std::vector<int>>(800);
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

  if (fExportSPWires) {
    art::Handle<std::vector<recob::Wire>> wires;
    e.getByLabel(fSPLabel, wires);

    fWiresZ = std::vector<std::vector<float>>(480);
    fWiresU = std::vector<std::vector<float>>(800);
    fWiresV = std::vector<std::vector<float>>(800);
    for (const recob::Wire& wire : *wires) {
      if (fGeom->ChannelToROP(wire.Channel()) == rIDZ) {
        fWiresZ[wire.Channel() - firstChZ] = std::vector<float>(4492);
        fillWires(wire, fWiresZ[wire.Channel() - firstChZ]);
      }
      else if (fGeom->ChannelToROP(wire.Channel()) == rIDU) {
        fWiresU[wire.Channel() - firstChU] = std::vector<float>(4492);
        fillWires(wire, fWiresU[wire.Channel() - firstChU]);
      }
      else if (fGeom->ChannelToROP(wire.Channel()) == rIDV) {
        fWiresV[wire.Channel() - firstChV] = std::vector<float>(4492);
        fillWires(wire, fWiresV[wire.Channel() - firstChV]);
      }
    }
  }

  if (fExportHits) {
    art::Handle<std::vector<recob::Hit>> hits;
    e.getByLabel(art::InputTag("gaushit"), hits);

    for (const recob::Hit& hit : *hits) {
      if (fGeom->ChannelToROP(hit.Channel()) == rIDU) {
        fHitROPs.push_back(0);
        fHitChs.push_back((int)(hit.Channel() - firstChU));
      }
      else if (fGeom->ChannelToROP(hit.Channel()) == rIDV) {
        fHitROPs.push_back(1);
        fHitChs.push_back((int)(hit.Channel() - firstChV));
      }
      else if (fGeom->ChannelToROP(hit.Channel()) == rIDZ) {
        fHitROPs.push_back(2);
        fHitChs.push_back((int)(hit.Channel() - firstChZ));
      }
      fHitPeakTimes.push_back(hit.PeakTime());
      fHitStartTicks.push_back(hit.StartTick());
      fHitEndTicks.push_back(hit.EndTick());
      fHitPeakSigmas.push_back(hit.SigmaPeakTime());
      fHitRMSs.push_back(hit.RMS());
      fHitPeakAmplitudes.push_back(hit.PeakAmplitude());
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
  fWiresZ.clear();
  fWiresU.clear();
  fWiresV.clear();
  fHitROPs.clear();
  fHitChs.clear();
  fHitPeakTimes.clear();
  fHitStartTicks.clear();
  fHitEndTicks.clear();
  fHitPeakSigmas.clear();
  fHitRMSs.clear();
  fHitPeakAmplitudes.clear();
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

void extrapolation::ExportFDSimRecoData::fillWires(
  const recob::Wire &wire, std::vector<float> &wireVec
)
{
  const std::vector<float> sigVec = wire.Signal();
  for (unsigned int tick = 0; tick < wireVec.size(); tick++) {
    wireVec[tick] = (float)sigVec[tick];
  }
}

DEFINE_ART_MODULE(extrapolation::ExportFDSimRecoData)
