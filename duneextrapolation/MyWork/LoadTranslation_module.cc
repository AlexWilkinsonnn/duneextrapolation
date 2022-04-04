////////////////////////////////////////////////////////////////////////
// Class:       LoadTranslation
// Plugin Type: producer (Unknown Unknown)
// File:        LoadTranslation_module.cc
//
// Generated Mon Mar 28 2022 by Alexander Wilkinson
//
// Loads true rawdigits, predicted rawdigits from model that takes ND
// response as input, and hits made from ND response.
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
  class LoadTranslation;

  struct digs;
  struct packet;
}

struct extrapolation::digs {
  Int_t ch;
  std::vector<short> digvec; // use Short_t, will need to change this in export_to_root.py also
};

struct extrapolation::packet {
  Int_t ch;
  Int_t tick;
  Int_t adc;
};

class extrapolation::LoadTranslation : public art::EDProducer {
public:
  explicit LoadTranslation(fhicl::ParameterSet const& p);

  LoadTranslation(LoadTranslation const&) = delete;
  LoadTranslation(LoadTranslation&&) = delete;
  LoadTranslation& operator=(LoadTranslation const&) = delete;
  LoadTranslation& operator=(LoadTranslation&&) = delete;

  void produce(art::Event& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  unsigned int fCIndex;
  unsigned int fTIndex;
  unsigned int fPIndex;
  
  int                            fEntry;
  int                            fNEntries;
  TTree*                         fTree;
  std::vector<int>*              fChannels;
  std::vector<std::vector<int>>* fTranslatedDigs;
  std::vector<std::vector<int>>* fTrueDigs;
  std::vector<std::vector<int>>* fNDPacket;

  readout::ROPID fRID;
  std::string fInputFileLoc;
};

extrapolation::LoadTranslation::LoadTranslation(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fCIndex (p.get<unsigned int>("CryoIndex")),
    fTIndex (p.get<unsigned int>("TpcIndex")),
    fPIndex (p.get<unsigned int>("PlaneIndex")),
    fInputFileLoc (p.get<std::string>("InputFileLoc"))
{
  produces<std::vector<raw::RawDigit>>("NetworkTranslated");
  produces<std::vector<raw::RawDigit>>("TrueTranslated");
  produces<std::vector<recob::Hit>>("NDPackets");
}

void extrapolation::LoadTranslation::produce(art::Event& e)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
  auto digsTranslated = std::make_unique<std::vector<raw::RawDigit>>();
  auto digsTrue = std::make_unique<std::vector<raw::RawDigit>>();
  auto hitsND = std::make_unique<std::vector<recob::Hit>>();

  if (fEntry < fNEntries) {
    fTree->GetEntry(fEntry);

    for (unsigned int ch = 0; ch < fGeom->Nchannels(); ++ch) {
      auto iCh = std::find(fChannels->begin(), fChannels->end(), (int)ch);

      if (iCh != fChannels->end()) {
        int chLocal = iCh - fChannels->begin();

        raw::RawDigit::ADCvector_t adcVec(4492);
        for (int tick = 0; tick < 4492; tick++){
          adcVec[tick] = (short)(*fTranslatedDigs)[chLocal][tick] + (short)900; // doing default collection pedestal manually for now just because I am afraid WC will have something hardcoded
        }
        raw::RawDigit rawDig(ch, adcVec.size(), adcVec);
        rawDig.SetPedestal(900);
        digsTranslated->push_back(rawDig);

        adcVec = raw::RawDigit::ADCvector_t(4492);
        for (int tick = 0; tick < 4492; tick++){
          adcVec[tick] = (short)(*fTrueDigs)[chLocal][tick] + (short)900;
        }
        rawDig = raw::RawDigit(ch, adcVec.size(), adcVec);
        rawDig.SetPedestal(900);
        digsTrue->push_back(rawDig); 
      }
      else {
        raw::RawDigit::ADCvector_t adcVec(4492, 900);
        raw::RawDigit rawDig(ch, adcVec.size(), adcVec);
        rawDig.SetPedestal(900);
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
  }

  fEntry++;
}

void extrapolation::LoadTranslation::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  const geo::CryostatID cID(fCIndex);
  const geo::TPCID tID(cID, fTIndex);
  const geo::PlaneID pID(tID, fPIndex);
  std::cout << pID << "\n";
  fRID = fGeom->WirePlaneToROP(pID);

  std::cout << "Reading from file " << fInputFileLoc << "\n";
  fChannels = nullptr;
  fTranslatedDigs = nullptr;
  fTrueDigs = nullptr;
  fNDPacket = nullptr;
  TFile* f = new TFile(fInputFileLoc.c_str());
  fTree = (TTree*)f->Get("digs_hits");
  fTree->SetBranchAddress("channels", &fChannels);
  fTree->SetBranchAddress("rawdigits_translated", &fTranslatedDigs);
  fTree->SetBranchAddress("rawdigits_true", &fTrueDigs);
  fTree->SetBranchAddress("nd_packets", &fNDPacket);

  fEntry = 0;
  fNEntries = fTree->GetEntries();
}

void extrapolation::LoadTranslation::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::LoadTranslation)
