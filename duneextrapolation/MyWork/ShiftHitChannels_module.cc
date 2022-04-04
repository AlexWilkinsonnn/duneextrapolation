////////////////////////////////////////////////////////////////////////
// Class:       ShiftHitChannels
// Plugin Type: producer (Unknown Unknown)
// File:        ShiftHitChannels_module.cc
//
// Generated Mon Mar 28 2022 by Alexander Wilkinson
//
// Shift ndpacket hits back to where they should be (fixing mistake I made). 
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TInterpreter.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace extrapolation {
  class ShiftHitChannels;
}

class extrapolation::ShiftHitChannels : public art::EDProducer {
public:
  explicit ShiftHitChannels(fhicl::ParameterSet const& p);

  ShiftHitChannels(ShiftHitChannels const&) = delete;
  ShiftHitChannels(ShiftHitChannels&&) = delete;
  ShiftHitChannels& operator=(ShiftHitChannels const&) = delete;
  ShiftHitChannels& operator=(ShiftHitChannels&&) = delete;

  void produce(art::Event& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  std::string fNDPacketLabel;
  int fChannelShift;
};

extrapolation::ShiftHitChannels::ShiftHitChannels(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fNDPacketLabel (p.get<std::string>("NDPacketLabel")),
  fChannelShift  (p.get<int>("ChannelShift"))
{
  produces<std::vector<recob::Hit>>("NDPacketsFixed");
}

void extrapolation::ShiftHitChannels::produce(art::Event& e)
{
  auto hitsNDFixed = std::make_unique<std::vector<recob::Hit>>();
  const auto hitsND = e.getValidHandle<std::vector<recob::Hit>> (fNDPacketLabel);

  for (const recob::Hit hit : *hitsND) {
    raw::ChannelID_t channel = (int)hit.Channel() + fChannelShift;
    geo::View_t view = fGeom->View(fGeom->ChannelToROP(channel));
    geo::WireID wireID = geo::WireID();
    for (geo::WireID wire : fGeom->ChannelToWire(channel)) {
      if (fGeom->View(wire.parentID()) == view) {
        wireID = wire;
        break;
      }
    }

    recob::Hit fixedHit(channel, hit.StartTick(), hit.EndTick(), hit.PeakTime(),
      hit.SigmaPeakTime(), hit.RMS(), hit.PeakAmplitude(), hit.SigmaPeakAmplitude(), hit.SummedADC(),
      hit.Integral(), hit.SigmaIntegral(), hit.Multiplicity(), hit.LocalIndex(), hit.GoodnessOfFit(),
      hit.DegreesOfFreedom(), view, hit.SignalType(), wireID); 
    
    hitsNDFixed->push_back(fixedHit);
  }

  e.put(std::move(hitsNDFixed), "NDPacketsFixed");
}

void extrapolation::ShiftHitChannels::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
  
  std::cout << "Shifting channel by " << fChannelShift << " for hits in " << fNDPacketLabel << "\n";
}

void extrapolation::ShiftHitChannels::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::ShiftHitChannels)
