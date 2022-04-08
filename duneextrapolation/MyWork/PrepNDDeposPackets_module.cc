////////////////////////////////////////////////////////////////////////
// Class:       PrepNDDeposPackets
// Plugin Type: producer (Unknown Unknown)
// File:        PrepNDDeposPackets_module.cc
//
// Generated at Thu Apr 07 2022 by Alexander Wilkinson.
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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace extrapolation {
  class PrepNDDeposPackets;
}

class extrapolation::PrepNDDeposPackets : public art::EDProducer {
public:
  explicit PrepNDDeposPackets(fhicl::ParameterSet const& p);

  PrepNDDeposPackets(PrepNDDeposPackets const&) = delete;
  PrepNDDeposPackets(PrepNDDeposPackets&&) = delete;
  PrepNDDeposPackets& operator=(PrepNDDeposPackets const&) = delete;
  PrepNDDeposPackets& operator=(PrepNDDeposPackets&&) = delete;

  void produce(art::Event& e) override;

  void beginJob() override;
  void endJob() override;

  void reset();
private:
  const geo::GeometryCore* fGeom;

  int                               fNEntries;
  int                               fEntry;
  TTree*                            fTreeDeposPackets;
  std::vector<std::vector<double>>* fDepos;
  std::vector<std::vector<double>>* fPackets;
  std::vector<double>*              fVertex;

  TTree*                           fTreePacketProjections;
  std::vector<std::vector<double>> fPacketProjection;

  std::string fNDDataLoc;
  unsigned int fCIndex;
  unsigned int fTIndex;
  int fEventNumber;
  double fYShift;
  double fZShift;
  double fTickShiftZ;
  double fTickShiftU;
  double fTickShiftV;
  geo::TPCID   fTID;
  geo::PlaneID fPIDZ;
  geo::PlaneID fPIDU;
  geo::PlaneID fPIDV;
};

extrapolation::PrepNDDeposPackets::PrepNDDeposPackets(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fNDDataLoc  (p.get<std::string>("NDDataLoc")),
    fCIndex     (p.get<unsigned int>("CryoIndex")),
    fTIndex     (p.get<unsigned int>("TpcIndex")),
    fYShift     (p.get<double>("YShift")),
    fZShift     (p.get<double>("ZShift")),
    fTickShiftZ (p.get<double>("TickShiftZ")),
    fTickShiftU (p.get<double>("TickShiftU")),
    fTickShiftV (p.get<double>("TickShiftV"))
{
  produces<std::vector<sim::SimEnergyDeposit>>();
  produces<std::vector<sim::SimEnergyDeposit>>("EventNumber");

  art::ServiceHandle<art::TFileService> tfs;

  fTreePacketProjections = tfs->make<TTree>("packet_projections", "packet_projections");
  fTreePacketProjections->Branch("projection", &fPacketProjection);
}

void extrapolation::PrepNDDeposPackets::produce(art::Event& e)
{
  this ->reset();

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
  // std::cout << detProp.ConvertTicksToX(2000, fPIDZ) << "\n"; = 163.705

  if (fEntry < fNEntries) { 
    fTreeDeposPackets->GetEntry(fEntry);

    // Get Vertex alignment info
    double vtxX = fVertex->at(2); // ND z is FD x (drift coordinate)
    double vtxTick = detProp.ConvertXToTicks(vtxX, fPIDZ); // Align wrt 2000 Z plane ticks
    double tickDiff = 2000 - vtxTick;
    double XShift = tickDiff * detProp.GetXTicksCoefficient(fTID.TPC, fTID.Cryostat);

    // Make and add the ND depos
    auto SEDs = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
    auto evNum = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

    for (const std::vector<double>& depo : *fDepos) {
      int trackID = (int)depo[0];
      int pdg = (int)depo[1];
      double zMin = depo[2] + fZShift;
      double zMax = depo[3] + fZShift;
      double yMin = depo[4] + fYShift;
      double yMax = depo[5] + fYShift;
      double xMin = depo[6] + XShift;
      double xMax = depo[7] + XShift;
      double tMin = depo[8];
      double tMax = depo[9];
      int electrons = (int)depo[10];
      double dE = depo[11];

      geo::Point_t posStart = geo::Point_t(xMin, yMin, zMin);
      geo::Point_t posEnd = geo::Point_t(xMax, yMax, zMax);

      sim::SimEnergyDeposit SED = sim::SimEnergyDeposit(
          0, electrons, 0, dE, posStart, posEnd, tMin, tMax, trackID, pdg); 
      SEDs->push_back(SED);

      // SED that stores an ID for this event
      sim::SimEnergyDeposit ID = sim::SimEnergyDeposit(0,0,0,0,posStart,posEnd,0,0,fEventNumber);
      evNum->push_back(ID);
    }

    e.put(std::move(SEDs));
    e.put(std::move(evNum), "EventNumber");

    // Get the wires and ticks of the ND packets using the vertex alignment
    for (const std::vector<double>& packet : *fPackets) {
      double z = packet[0] + fZShift;
      double y = packet[1] + fYShift;
      double x = packet[2] + XShift;

      geo::Point_t packetLoc(x, y, z);

      raw::ChannelID_t chZ = fGeom->NearestChannel(packetLoc, fPIDZ);
      raw::ChannelID_t chU = fGeom->NearestChannel(packetLoc, fPIDU);
      raw::ChannelID_t chV = fGeom->NearestChannel(packetLoc, fPIDV);

      double tickRawZ = detProp.ConvertXToTicks(x, fPIDZ);
      tickRawZ -= fTickShiftZ;
      unsigned int tickZ = (unsigned int)tickRawZ;
      double tickRawU = detProp.ConvertXToTicks(x, fPIDU);
      tickRawU -= fTickShiftU;
      unsigned int tickU = (unsigned int)tickRawZ;
      double tickRawV = detProp.ConvertXToTicks(x, fPIDV);
      tickRawV -= fTickShiftV;
      unsigned int tickV = (unsigned int)tickRawZ;

      std::vector<double> projection(10, 0.0);
      projection[0] = fEventNumber;
      projection[1] = z;
      projection[2] = y;
      projection[3] = x;
      projection[4] = chZ;
      projection[5] = tickZ;
      projection[6] = chU;
      projection[7] = tickU;
      projection[8] = chV;
      projection[9] = tickV;
      fPacketProjection.push_back(projection);
    }

    fTreePacketProjections->Fill();

    fEventNumber++;
  }
  else {
    std::cout << "Gone beyond number of entries in tree (" << fNEntries << ")\n";
    fEntry++;
  }
}
void extrapolation::PrepNDDeposPackets::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  const geo::CryostatID cID(fCIndex);
  fTID = geo::TPCID(cID, fTIndex);
  for (geo::PlaneID pID : fGeom->IteratePlaneIDs(fTID)) {
    if (fGeom->View(pID) == geo::kZ) {
      std::cout << "Z plane: " << pID << "\n";
      fPIDZ = pID;
    }
    else if (fGeom->View(pID) == geo::kU) {
      std::cout << "U plane: " << pID << "\n";
      fPIDU = pID;
    }
    else if (fGeom->View(pID) == geo::kV) {
      std::cout << "V plane: " << pID << "\n";
      fPIDV = pID;
    }
  }

//   fYShift = 200; // cm
//   fZShift = -14.8615; // cm 

  fDepos = nullptr;
  fPackets = nullptr;
  fVertex = nullptr;
  TFile* f = new TFile(fNDDataLoc.c_str());
  fTreeDeposPackets = (TTree*)f->Get("ND_depos_packets");
  fTreeDeposPackets->SetBranchAddress("nd_depos", &fDepos);
  fTreeDeposPackets->SetBranchAddress("nd_packets", &fPackets);
  fTreeDeposPackets->SetBranchAddress("vertex", &fVertex);

  fEventNumber = 0;
  fNEntries = fTreeDeposPackets->GetEntries();
}

void extrapolation::PrepNDDeposPackets::endJob()
{
}

void extrapolation::PrepNDDeposPackets::reset()
{
  fPacketProjection.clear();
}

DEFINE_ART_MODULE(extrapolation::PrepNDDeposPackets)

