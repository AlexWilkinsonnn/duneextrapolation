////////////////////////////////////////////////////////////////////////
// Class:       CheckBackTrackingSCWires
// Plugin Type: analyzer (Unknown Unknown)
// File:        CheckBackTrackingSCWires_module.cc
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
#include "lardataobj/RecoBase/Wire.h"

#include <map>

namespace extrapolation {
  class CheckBackTrackingSCWires;
}

class extrapolation::CheckBackTrackingSCWires : public art::EDAnalyzer {
public:
  explicit CheckBackTrackingSCWires(fhicl::ParameterSet const& p);

  CheckBackTrackingSCWires(CheckBackTrackingSCWires const&) = delete;
  CheckBackTrackingSCWires(CheckBackTrackingSCWires&&) = delete;
  CheckBackTrackingSCWires& operator=(CheckBackTrackingSCWires const&) = delete;
  CheckBackTrackingSCWires& operator=(CheckBackTrackingSCWires&&) = delete;

  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

  void reset();
private:
  const geo::GeometryCore* fGeom;

  bool         fSpecifyCh;
  unsigned int fTPCIndex;
  unsigned int fPlaneIndex;
  unsigned int fChLocal;

  std::string fSCLabel;
  std::string fWireLabel;

  unsigned int fSpecificChGlobal;

  TTree*                           fTreeWires;
  std::vector<std::vector<double>> fWires;
  std::vector<int>                 fWiresChs;
  std::vector<std::vector<double>> fSC;
  std::vector<int>                 fSCChs;
  std::map<int, int>               fChTypes; // 0 induction, 1 collection
  std::vector<double>              fSpecificWire;
  std::vector<double>              fSpecificSC;
  std::vector<int>                 fSpecificTPCPlaneCh;
};


extrapolation::CheckBackTrackingSCWires::CheckBackTrackingSCWires(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSpecifyCh  (p.get<bool>         ("SpecifyChannel")),
    fTPCIndex   (p.get<unsigned int> ("TPCIndex")),
    fPlaneIndex (p.get<unsigned int> ("PlaneIndex")),
    fChLocal    (p.get<unsigned int> ("ChannelLocal")),
    fSCLabel    (p.get<std::string>  ("SCLabel")),
    fWireLabel  (p.get<std::string>  ("WireLabel"))
{
  consumes<std::vector<sim::SimChannel>>(fSCLabel);
  consumes<std::vector<recob::Wire>>(fWireLabel);

  art::ServiceHandle<art::TFileService> tfs;

  fTreeWires = tfs->make<TTree>("sc_wire", "sc_wire");
  fTreeWires->Branch("simchannels", &fSC);
  fTreeWires->Branch("simchannels_chs", &fSCChs);
  fTreeWires->Branch("wires", &fWires);
  fTreeWires->Branch("wires_chs", &fWiresChs);
  fTreeWires->Branch("ch_types", &fChTypes);
  if (fSpecifyCh) {
    fTreeWires->Branch("specific_wire", &fSpecificWire);
    fTreeWires->Branch("specific_simchannel", &fSpecificSC);
    fTreeWires->Branch("specific_tpc_plane_ch", &fSpecificTPCPlaneCh);
  }
}

void extrapolation::CheckBackTrackingSCWires::analyze(art::Event const& e)
{
  this->reset();

  art::Handle<std::vector<sim::SimChannel>> simChs;
  e.getByLabel(fSCLabel, simChs);
  art::Handle<std::vector<recob::Wire>> wires;
  e.getByLabel(fWireLabel, wires);

  for (sim::SimChannel simCh : *simChs) {
    raw::ChannelID_t ch = simCh.Channel();
    std::vector<double> respVec(6000, 0.0);
    for (sim::TDCIDE tickIde : simCh.TDCIDEMap()) {
      int tick = (int)tickIde.first;
      for (sim::IDE ide : tickIde.second) {
        respVec[tick] += (double)ide.energy;
      }
    }
    if (std::accumulate(respVec.begin(), respVec.end(), 0.0) != 0.0) {
      fChTypes[(int)ch] = fGeom->SignalType(ch);
      fSC.push_back(respVec);
      fSCChs.push_back((int)ch);
    }
  }

  for (recob::Wire wire : *wires) {
    raw::ChannelID_t ch = wire.Channel();
    std::vector<double> respVec(6000, 0.0);
    int tick = 0;
    for (float mag : wire.Signal()) {
      if (mag != 0) {
        respVec[tick] += (double)mag;
      }
      tick++;
    }
    if (std::accumulate(respVec.begin(), respVec.end(), 0.0) != 0.0) {
      fChTypes[(int)ch] = fGeom->SignalType(ch);
      fWires.push_back(respVec);
      fWiresChs.push_back((int)ch);
    }
  }

  fTreeWires->Fill();
}

void extrapolation::CheckBackTrackingSCWires::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  // // XXX In progess
  // if (fSpecifyCh) {
  //   const geo::CryostatID cID(0);
  //   const geo::TPCID tID(cID, fTPCIndex);
  //   const geo::PlaneID pID(tID, fPlaneIndex);
  //   const readout::ROPID rID = fGeom->WirePlaneToROP(pID);
  //   fSpecificChGlobal = fChLocal + fGeom->FirstChannelInROP(rID);
  //   std::cout << "TPC " << fTPCIndex
  //             << " Plane " << fPlaneIndex << " (" <<
  //             fGeom->SignalType(pID) == geo::SigType_t::kCollection ? "collection" : "induction"
  //             << ") ch " << fChLocal << " is global ch " << fSpecificChGlobal << "\n";
  // }
}

void extrapolation::CheckBackTrackingSCWires::endJob()
{
}

void extrapolation::CheckBackTrackingSCWires::reset()
{
  fSC.clear();
  fSCChs.clear();
  fWires.clear();
  fWiresChs.clear();
  fChTypes.clear();
  fSpecificWire.clear();
  fSpecificSC.clear();
  fSpecificTPCPlaneCh.clear();
}

DEFINE_ART_MODULE(extrapolation::CheckBackTrackingSCWires)

