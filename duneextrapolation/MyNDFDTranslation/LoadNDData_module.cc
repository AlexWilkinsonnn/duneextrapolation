////////////////////////////////////////////////////////////////////////
// Class:       LoadNDData
// Plugin Type: producer (Unknown Unknown)
// File:        LoadNDData_module.cc
//
// Copied on Mon 20 Jun 22 by Alex Wilkinson.
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
#include "lardataobj/RawData/raw.h"

namespace extrapolation {
  class LoadNDData;
}

class extrapolation::LoadNDData : public art::EDProducer {
public:
  explicit LoadNDData(fhicl::ParameterSet const& p);

  LoadNDData(LoadNDData const&) = delete;
  LoadNDData(LoadNDData&&) = delete;
  LoadNDData& operator=(LoadNDData const&) = delete;
  LoadNDData& operator=(LoadNDData&&) = delete;

  void produce(art::Event& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  // Reading input ND tree
  int                               fNEntries;
  int                               fEntry;
  TTree*                            fTreeDeposPackets;
  std::vector<std::vector<double>>* fDepos;
  std::vector<std::vector<double>>* fPackets;
  // NOTE Need to decide how to keep track of which events are the same between larsoft and the ND
  // reco chain. An int wont be good enough but but fine for now. Currently am just goint to use
  // the entry number in the tree but will eventually need to have it created in the ND chain and
  // read into here throught the ND input tree.
  // int                               fNDEventID;

  // fhicl params
  std::string  fNDDataLoc;
  unsigned int fCIndex;
  unsigned int fTIndex;
  double       fXShift;
  double       fYShift;
  double       fZShift;
  bool         fLoadTrueDepos;
};

extrapolation::LoadNDData::LoadNDData(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fNDDataLoc            (p.get<std::string>("NDDataLoc")),
    fXShift               (p.get<double>("XShift")),
    fYShift               (p.get<double>("YShift")),
    fZShift               (p.get<double>("ZShift")),
    fLoadTrueDepos        (p.get<bool>("LoadTrueDepos"))
{
  produces<std::vector<sim::SimEnergyDeposit>>("NDPackets");
  produces<std::vector<sim::SimEnergyDeposit>>("EventID");

  if (fLoadTrueDepos) {
    produces<std::vector<sim::SimEnergyDeposit>>("NDDepos");
  }
}

void extrapolation::LoadNDData::produce(art::Event& e)
{
  if (fEntry < fNEntries) {
    fTreeDeposPackets->GetEntry(fEntry);

    // Make and add the ND depos
    if (fLoadTrueDepos) {
      auto NDDepos = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

      for (const std::vector<double>& depo : *fDepos) {
        int trackID = (int)depo[0];
        int pdg = (int)depo[1];
        double zMin = depo[2] + fZShift;
        double zMax = depo[3] + fZShift;
        double yMin = depo[4] + fYShift;
        double yMax = depo[5] + fYShift;
        double xMin = depo[6] + fXShift;
        double xMax = depo[7] + fXShift;
        double tMin = depo[8];
        double tMax = depo[9];
        int electrons = (int)depo[10];
        double dE = depo[11];

        geo::Point_t posStart = geo::Point_t(xMin, yMin, zMin);
        geo::Point_t posEnd = geo::Point_t(xMax, yMax, zMax);

        sim::SimEnergyDeposit NDDepo = sim::SimEnergyDeposit(
            0, electrons, 0, dE, posStart, posEnd, tMin, tMax, trackID, pdg);
        NDDepos->push_back(NDDepo);
      }

      e.put(std::move(NDDepos), "NDDepos");
    }

    // SED that stores an ID for this event.
    // TODO Implement a proper ND event ID.
    auto evID = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

    geo::Point_t posStart = geo::Point_t(0,0,0);
    geo::Point_t posEnd = geo::Point_t(0,0,0);
    sim::SimEnergyDeposit ID = sim::SimEnergyDeposit(0,0,0,0, posStart, posEnd, 0,0, fEntry);
    evID->push_back(ID);

    e.put(std::move(evID), "EventID");

    // Write NDPackets into event. Leave projection to the pixelmap maker.
    auto NDPackets = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

    for (const std::vector<double>& packet : *fPackets) {
      double z = packet[0] + fZShift;
      double y = packet[1] + fYShift;
      double x = packet[2] + fXShift;
      double adc = packet[4];
      double NDDrift = packet[5];

      geo::Point_t posStart(x, y, z);
      geo::Point_t posEnd(x, y, z);

      // Hijacking SimEnergyDeposit to store packets for now.
      // electrons <-> adc, scintYieldRatio <-> NDDrift.
      // TODO Define my own data product to store these.
      sim::SimEnergyDeposit NDPacket = sim::SimEnergyDeposit(
        0, int(adc), 0, NDDrift, posStart, posEnd, 0,0,0,0);
      NDPackets->push_back(NDPacket);
    }

    e.put(std::move(NDPackets), "NDPackets");
  }
  else {
    std::cout << "Gone beyond number of entries in tree (" << fNEntries << ")\n"
      << "Empty events are being written out\n";
  }

  fEntry++;
}

void extrapolation::LoadNDData::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  std::cout << "Reading file from " << fNDDataLoc << "\n";
  fEntry = 0;
  fDepos = nullptr;
  fPackets = nullptr;
  TFile* f = new TFile(fNDDataLoc.c_str());
  fTreeDeposPackets = (TTree*)f->Get("ND_depos_packets");
  fTreeDeposPackets->SetBranchAddress("nd_depos", &fDepos);
  fTreeDeposPackets->SetBranchAddress("nd_packets", &fPackets);

  fNEntries = fTreeDeposPackets->GetEntries();
  std::cout << "File has " << fNEntries << " entries\n";
}

void extrapolation::LoadNDData::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::LoadNDData)

