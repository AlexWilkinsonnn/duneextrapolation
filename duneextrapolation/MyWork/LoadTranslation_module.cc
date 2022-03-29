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
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

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
  
  int                                     fEntry;
  int                                     fNEntries;
  TTree*                                  fTree;
  std::vector<int>                        fChannels;
  std::vector<std::vector<int>>           fTranslatedDigs;
  std::vector<std::vector<int>>           fTrueDigs;
  std::vector<std::map<std::string, int>> fNDPacket;

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
  produces<std::vector<raw::RawDigit>>("Translated");
  produces<std::vector<raw::RawDigit>>("True");
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
      auto iCh = std::find(fChannels.begin(), fChannels.end(), (int)ch);

      if (iCh != fChannels.end()) {
        int chLocal = iCh - fChannels.begin();

        raw::RawDigit::ADCvector_t adcVec(4492);
        for (int tick = 0; tick < 4492; tick++){
          adcVec[tick] = (short)fTranslatedDigs[chLocal][tick];
        }
        raw::RawDigit rawDig(ch, adcVec.size(), adcVec);
        rawDig.SetPedestal(0);
        digsTranslated->push_back(rawDig);

        adcVec = raw::RawDigit::ADCvector_t(4492);
        for (int tick = 0; tick < 4492; tick++){
          adcVec[tick] = (short)fTrueDigs[chLocal][tick];
        }
        rawDig = raw::RawDigit(ch, adcVec.size(), adcVec);
        rawDig.SetPedestal(0);
        digsTrue->push_back(rawDig); 
      }
      else {
        raw::RawDigit::ADCvector_t adcVec(4492, 0);
        raw::RawDigit rawDig(ch, adcVec.size(), adcVec);
        rawDig.SetPedestal(0);
        digsTrue->push_back(rawDig);
        digsTranslated->push_back(rawDig);   
      }
    }

    for (std::map<std::string, int> ndPacket : fNDPacket) {
      float peakTime = ndPacket["tick"];
      float integral = ndPacket["adc"];
      raw::ChannelID_t channel = ndPacket["ch"];
      geo::View_t view = fGeom->View(fGeom->ChannelToROP(ndPacket["ch"]));

      geo::WireID wireID = geo::WireID();
      for (geo::WireID wire : fGeom->ChannelToWire(ndPacket["ch"])) { 
        if (fGeom->View(wire.parentID()) == view) {
          wireID = wire;
          break;
        }
      }

      recob::Hit hit(channel, raw::TDCtick_t(0), raw::TDCtick_t(0), peakTime, float(-1.0),
        float(0.0), float(0.0), float(-1.0), float(0.0), integral, float(-1.0), short(0),
        short(-1), float(0.0), int(-1), view, geo::SigType_t(geo::kMysteryType), wireID);

      hitsND->push_back(hit);
    }

    e.put(std::move(digsTrue), "True");
    e.put(std::move(digsTranslated), "Translated");
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

  TFile* f = new TFile(fInputFileLoc.c_str());
  fTree = (TTree*)f->Get("tree");
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
