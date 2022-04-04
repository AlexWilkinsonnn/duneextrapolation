////////////////////////////////////////////////////////////////////////
// Class:       ScaleNDHits
// Plugin Type: producer (Unknown Unknown)
// File:        ScaleNDHits_module.cc
//
// Generated Mon Apr 4 2022 by Alexander Wilkinson
//
// Scale ND hits using the means of the FD and ND hits.
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
  class ScaleNDHits;
}

class extrapolation::ScaleNDHits : public art::EDProducer {
public:
  explicit ScaleNDHits(fhicl::ParameterSet const& p);

  ScaleNDHits(ScaleNDHits const&) = delete;
  ScaleNDHits(ScaleNDHits&&) = delete;
  ScaleNDHits& operator=(ScaleNDHits const&) = delete;
  ScaleNDHits& operator=(ScaleNDHits&&) = delete;

  void produce(art::Event& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  std::string fNDPacketLabel;
  double fNDHitsMean;
  double fFDTrueHitsMean;
};

extrapolation::ScaleNDHits::ScaleNDHits(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fNDPacketLabel  (p.get<std::string>("NDPacketLabel")),
  fNDHitsMean     (p.get<double>("FDTrueHitsMean")),
  fFDTrueHitsMean (p.get<double>("FDTrueHitsMean"))
{
  produces<std::vector<recob::Hit>>("NDPacketsScaled");
  
  consumes<std::vector<recob::Hit>>(fNDPacketLabel);
}

void extrapolation::ScaleNDHits::produce(art::Event& e)
{
  auto hitsNDScaled = std::make_unique<std::vector<recob::Hit>>();
  const auto hitsND = e.getValidHandle<std::vector<recob::Hit>> (fNDPacketLabel);

  for (const recob::Hit hit : *hitsND) {
    raw::ChannelID_t channel = (int)hit.Channel();
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
    
    hitsNDScaled->push_back(fixedHit);
  }

  e.put(std::move(hitsNDScaled), "NDPacketsScaled");
}

void extrapolation::ScaleNDHits::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
  
  // std::cout << "Shifting channel by " << fChannelShift << " for hits in " << fNDPacketLabel << "\n";
}

void extrapolation::ScaleNDHits::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::ScaleNDHits)
