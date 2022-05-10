////////////////////////////////////////////////////////////////////////
// Class:       LoadTranlsationKeepSED
// Plugin Type: filter (Unknown Unknown)
// File:        LoadTranlsationKeepSED_module.cc
//
// Copied May 10 2022 by Alexander Wilkinson
//
// Load tranalsation results into event if event ID has results.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
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
#include <map>
#include <set>

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
  class LoadTranlsationKeepSED;
}

class extrapolation::LoadTranlsationKeepSED : public art::EDFilter {
public:
  explicit LoadTranlsationKeepSED(fhicl::ParameterSet const& p);

  LoadTranlsationKeepSED(LoadTranlsationKeepSED const&) = delete;
  LoadTranlsationKeepSED(LoadTranlsationKeepSED&&) = delete;
  LoadTranlsationKeepSED& operator=(LoadTranlsationKeepSED const&) = delete;
  LoadTranlsationKeepSED& operator=(LoadTranlsationKeepSED&&) = delete;

  bool filter(art::Event& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  std::map<int, int> evNumToTreeEntry;

  TTree*                         fTree;
  std::vector<int>*              fChannels;
  std::vector<std::vector<int>>* fTranslatedDigs;
  std::vector<std::vector<int>>* fTrueDigs;
  std::vector<std::vector<int>>* fNDPacket;

  std::string fInputFileLoc;
  std::string fEvNumLabel;
  double      fNDHitsScaleFactor;
};

extrapolation::LoadTranlsationKeepSED::LoadTranlsationKeepSED(fhicl::ParameterSet const& p)
  : EDFilter{p},
    fInputFileLoc      (p.get<std::string>("InputFileLoc")),
    fEvNumLabel        (p.get<std::string>("EvNumLabel")),
    fNDHitsScaleFactor (p.get<double>("NDHitsScaleFactor")) 
{
  produces<std::vector<raw::RawDigit>>("NetworkTranslated");
  produces<std::vector<raw::RawDigit>>("TrueTranslated");
  produces<std::vector<recob::Hit>>("NDPackets");

  consumes<std::vector<sim::SimEnergyDeposit>>(fEvNumLabel);

  TFile* f = new TFile(fInputFileLoc.c_str());
  fTree = (TTree*)f->Get("digs_hits");
  const int nEntries = fTree->GetEntries();

  int evNum = -1;
  fTree->SetBranchAddress("ev_num", &evNum);

  for (int i = 0; i < nEntries; i++) { 
    evNum = -1;
    fTree->GetEntry(i); 

    std::cout << evNum << "\n";
    evNumToTreeEntry[evNum] = i;
  }
}

bool extrapolation::LoadTranlsationKeepSED::filter(art::Event& e)
{
  const auto evNumVec = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEvNumLabel);

  int evNum = evNumVec->at(0).TrackID();

  if (!evNumToTreeEntry.count(evNum)) {
    return false;
  }

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  auto digsTranslated = std::make_unique<std::vector<raw::RawDigit>>();
  auto digsTrue = std::make_unique<std::vector<raw::RawDigit>>();
  auto hitsND = std::make_unique<std::vector<recob::Hit>>();

  fTree->GetEntry(evNumToTreeEntry[evNum]);

  for (unsigned int ch = 0; ch < fGeom->Nchannels(); ++ch) {
    short ped = fGeom->View(fGeom->ChannelToROP(ch)) == geo::kZ ? 900 : 2350;

    auto iCh = std::find(fChannels->begin(), fChannels->end(), (int)ch);

    if (iCh != fChannels->end()) {
      int chLocal = iCh - fChannels->begin();

      raw::RawDigit::ADCvector_t adcVec(4492);
      for (int tick = 0; tick < 4492; tick++){
        adcVec[tick] = (short)(*fTranslatedDigs)[chLocal][tick] + ped;// doing default collection pedestal manually for now just because I am afraid WC will have something hardcoded
      }
      raw::RawDigit rawDig(ch, adcVec.size(), adcVec);
      rawDig.SetPedestal(ped);
      digsTranslated->push_back(rawDig);

      adcVec = raw::RawDigit::ADCvector_t(4492);
      for (int tick = 0; tick < 4492; tick++){
        adcVec[tick] = (short)(*fTrueDigs)[chLocal][tick] + ped;
      }
      rawDig = raw::RawDigit(ch, adcVec.size(), adcVec);
      rawDig.SetPedestal(ped);
      digsTrue->push_back(rawDig); 
    }
    else {
      raw::RawDigit::ADCvector_t adcVec(4492, ped);
      raw::RawDigit rawDig(ch, adcVec.size(), adcVec);
      rawDig.SetPedestal(ped);
      digsTrue->push_back(rawDig);
      digsTranslated->push_back(rawDig);   
    }
  }

  for (std::vector<int> ndPacket : *fNDPacket) {
    float peakTime = (float)ndPacket[1];
    float integral = (float)(ndPacket[2]*16);
    float summedADC = (float)(ndPacket[2]*16);
    float peakAmplitude = integral/5.0;
    float rms = 3.0; // fairly arbitrary choice
    raw::TDCtick_t startTick = (raw::TDCtick_t)(peakTime);
    raw::TDCtick_t endTick = (raw::TDCtick_t)(peakTime + 5.0);
    raw::ChannelID_t channel = (unsigned int)(ndPacket[0]);
    geo::View_t view = fGeom->View(fGeom->ChannelToROP((unsigned int)ndPacket[0]));

    // Empirical scaling to be ~ FD scale
    integral = (float)((double)integral * fNDHitsScaleFactor);

    geo::WireID wireID = geo::WireID();
    for (geo::WireID wire : fGeom->ChannelToWire((unsigned int)ndPacket[0])) {
      if (fGeom->View(wire.parentID()) == view) {
        wireID = wire;
        break;
      }
    }

    recob::Hit hit(channel, startTick, endTick, peakTime, float(-1.0),
      rms, peakAmplitude, float(-1.0), summedADC, integral, float(-1.0), short(0),
      short(-1), float(0.0), int(-1), view, geo::SigType_t(geo::kMysteryType), wireID);

    hitsND->push_back(hit);
  }

  e.put(std::move(digsTrue), "TrueTranslated");
  e.put(std::move(digsTranslated), "NetworkTranslated");
  e.put(std::move(hitsND), "NDPackets");

  return true;
}

void extrapolation::LoadTranlsationKeepSED::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
}

void extrapolation::LoadTranlsationKeepSED::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::LoadTranlsationKeepSED)
