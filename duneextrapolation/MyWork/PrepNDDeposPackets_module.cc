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

  // Reading input ND tree
  int                               fNEntries;
  int                               fEntry;
  TTree*                            fTreeDeposPackets;
  std::vector<std::vector<double>>* fDepos;
  std::vector<std::vector<double>>* fPackets;
  std::vector<double>*              fVertex;

  // Writing output ND projectiosn tree
  TTree*                           fTreePacketProjections;
  std::vector<std::vector<double>> fPacketProjection;
  int                              fEventID;
  std::vector<double>              fVertexOut;

  // fhicl params
  std::string  fNDDataLoc;
  unsigned int fCIndex;
  unsigned int fTIndex;
  double       fYShift;
  double       fZShift;
  double       fTickShiftZ;
  double       fTickShiftU;
  double       fTickShiftV;
  bool         fOnlyProjections;
  bool         fHighResProjection;

  // Other members
  int            fEventNumber;
  geo::TPCID     fTID;
  geo::PlaneID   fPIDZ;
  geo::PlaneID   fPIDU;
  geo::PlaneID   fPIDV;
  readout::ROPID fRIDZ;
  readout::ROPID fRIDU;
  readout::ROPID fRIDV;
};

extrapolation::PrepNDDeposPackets::PrepNDDeposPackets(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fNDDataLoc         (p.get<std::string>("NDDataLoc")),
    fCIndex            (p.get<unsigned int>("CryoIndex")),
    fTIndex            (p.get<unsigned int>("TpcIndex")),
    fYShift            (p.get<double>("YShift")),
    fZShift            (p.get<double>("ZShift")),
    fTickShiftZ        (p.get<double>("TickShiftZ")),
    fTickShiftU        (p.get<double>("TickShiftU")),
    fTickShiftV        (p.get<double>("TickShiftV")),
    fOnlyProjections   (p.get<bool>("OnlyProjections")),
    fHighResProjection (p.get<bool>("HighResProjection"))
{
  produces<std::vector<sim::SimEnergyDeposit>>();
  produces<std::vector<sim::SimEnergyDeposit>>("EventNumber");

  art::ServiceHandle<art::TFileService> tfs;

  fTreePacketProjections = tfs->make<TTree>("packet_projections", "packet_projections");
  fTreePacketProjections->Branch("projection", &fPacketProjection);
  fTreePacketProjections->Branch("eventid", &fEventID);
  fTreePacketProjections->Branch("vertex", &fVertexOut);
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

    if (!fOnlyProjections) {
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
      }
    }

    // SED that stores an ID for this event
    geo::Point_t posStart = geo::Point_t(0,0,0);
    geo::Point_t posEnd = geo::Point_t(0,0,0);
    sim::SimEnergyDeposit ID = sim::SimEnergyDeposit(0,0,0,0,posStart,posEnd,0,0,fEventNumber);
    evNum->push_back(ID);

    e.put(std::move(SEDs));
    e.put(std::move(evNum), "EventNumber");

    // Dont want to use this pair so dont write out the packet projections
    // Get the wires and ticks of the ND packets using the vertex alignment
    for (const std::vector<double>& packet : *fPackets) {
      double z = packet[0] + fZShift;
      double y = packet[1] + fYShift;
      double x = packet[2] + XShift;
      double adc = packet[4];
      double NDDrift = packet[5];

      geo::Point_t packetLoc(x, y, z);

      const geo::PlaneGeo pGeoZ = fGeom->Plane(fPIDZ);
      const geo::PlaneGeo pGeoU = fGeom->Plane(fPIDU);
      const geo::PlaneGeo pGeoV = fGeom->Plane(fPIDV);

      raw::ChannelID_t chZ = fGeom->NearestChannel(packetLoc, fPIDZ) - fGeom->FirstChannelInROP(fRIDZ);
      raw::ChannelID_t chU = fGeom->NearestChannel(packetLoc, fPIDU) - fGeom->FirstChannelInROP(fRIDU);
      raw::ChannelID_t chV = fGeom->NearestChannel(packetLoc, fPIDV) - fGeom->FirstChannelInROP(fRIDV);

      double tickRawZ = detProp.ConvertXToTicks(x, fPIDZ);
      tickRawZ += fTickShiftZ;
      unsigned int tickZ = (unsigned int)tickRawZ;
      double tickRawU = detProp.ConvertXToTicks(x, fPIDU);
      tickRawU += fTickShiftU;
      unsigned int tickU = (unsigned int)tickRawU;
      double tickRawV = detProp.ConvertXToTicks(x, fPIDV);
      tickRawV += fTickShiftV;
      unsigned int tickV = (unsigned int)tickRawV;

      double driftDistanceZ = pGeoZ.DistanceFromPlane(packetLoc);
      double driftDistanceU = pGeoU.DistanceFromPlane(packetLoc);
      double driftDistanceV = pGeoV.DistanceFromPlane(packetLoc);

      double wireCoordZ = pGeoZ.WireCoordinate(packetLoc);
      double wireDistanceZ = (wireCoordZ - (double)(int)(0.5 + wireCoordZ)) * pGeoZ.WirePitch();
      double wireCoordU = pGeoU.WireCoordinate(packetLoc);
      double wireDistanceU = (wireCoordU - (double)(int)(0.5 + wireCoordU)) * pGeoU.WirePitch();
      double wireCoordV = pGeoV.WireCoordinate(packetLoc);
      double wireDistanceV = (wireCoordV - (double)(int)(0.5 + wireCoordV)) * pGeoV.WirePitch();

      // Project to wires and ticks with a much higher resolution
      // Try 9 also
      if (fHighResProjection) {
        // std::cout << chZ << " -- ";
        const geo::WireID wIDZ = fGeom->NearestWireID(packetLoc, fPIDZ);
        // Cant have negative channel numbers, want the most negative pixel from the zeroth channel
        // to map to the zeroth pixel of the high res FD projection.
        const int wNoHighResShiftedZ = (int)(((wireCoordZ + 0.50) * 8.0));
        chZ = (raw::ChannelID_t)(wNoHighResShiftedZ - (wIDZ.Wire * 8) + (chZ * 8));
        if (chZ >= 3840 || chZ < 0) {
          std::cout << wNoHighResShiftedZ << " -- " << wIDZ.Wire * 8 << "\n";
          throw std::runtime_error("brokey");
        }
        // std::cout << chZ << " | ";
        // std::cout << tickZ << " -- ";
        tickZ = (unsigned int)((tickRawZ * 8.0) + 0.5);
        // std::cout << tickZ << " -- ";
        // std::cout << chU << " -- ";
        const geo::WireID wIDU = fGeom->NearestWireID(packetLoc, fPIDU);
        const int wNoHighResShiftedU = (int)(((wireCoordU + 0.50) * 8.0));
        chU = (raw::ChannelID_t)(wNoHighResShiftedU - (wIDU.Wire * 8) + (chU * 8));
        // My understanding is that this should not occur but it sometimes does so just going to fix it here
        if (chU  >= 6400 || chU < 0) {
          std::cout << wNoHighResShiftedU << " -- " << wIDU.Wire * 8 << "\n";
          throw std::runtime_error("brokey");
        }
        // std::cout << chU << " | ";
        // std::cout << tickU << " -- ";
        tickU = (unsigned int)((tickRawU * 8.0) + 0.5);
        // std::cout << tickU << "\n";
        // std::cout << chV << " -- ";
        const geo::WireID wIDV = fGeom->NearestWireID(packetLoc, fPIDV);
        const int wNoHighResShiftedV = (int)(((wireCoordV + 0.50) * 8.0));
        chV = (raw::ChannelID_t)(wNoHighResShiftedV - (wIDV.Wire * 8) + (chV * 8));
        if (chV >= 6400 || chV < 0) {
          std::cout << wNoHighResShiftedV << " -- " << wIDV.Wire * 8 << "\n";
          throw std::runtime_error("brokey");
        }
        // std::cout << chV << " | \n";
        // std::cout << tickV << " -- ";
        tickV = (unsigned int)((tickRawV * 8.0) + 0.5);
        // std::cout << tickV << "\n";
      }

      std::vector<double> projection(17, 0.0);
      projection[0] = z;
      projection[1] = y;
      projection[2] = x;
      projection[3] = chZ;
      projection[4] = tickZ;
      projection[5] = chU;
      projection[6] = tickU;
      projection[7] = chV;
      projection[8] = tickV;
      projection[9] = adc;
      projection[10] = NDDrift;
      projection[11] = driftDistanceZ;
      projection[12] = driftDistanceU;
      projection[13] = driftDistanceV;
      projection[14] = wireDistanceZ;
      projection[15] = wireDistanceU;
      projection[16] = wireDistanceV;
      fPacketProjection.push_back(projection);
    }

    // Store Event id number
    fEventID = fEventNumber;

    // Want to keep vertex info
    fVertexOut = *fVertex;

    fTreePacketProjections->Fill();

    fEventNumber++;
  }
  else {
    std::cout << "Gone beyond number of entries in tree (" << fNEntries << ")\n";
  }

  fEntry++;
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
      fRIDZ = fGeom->WirePlaneToROP(pID);
    }
    else if (fGeom->View(pID) == geo::kU) {
      std::cout << "U plane: " << pID << "\n";
      fPIDU = pID;
      fRIDU = fGeom->WirePlaneToROP(pID);
    }
    else if (fGeom->View(pID) == geo::kV) {
      std::cout << "V plane: " << pID << "\n";
      fPIDV = pID;
      fRIDV = fGeom->WirePlaneToROP(pID);
    }
  }

  std::cout << "Reading file from " << fNDDataLoc << "\n";
  fEntry = 0;
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
  std::cout << "File has " << fNEntries << " entries\n";
}

void extrapolation::PrepNDDeposPackets::endJob()
{
}

void extrapolation::PrepNDDeposPackets::reset()
{
  fPacketProjection.clear();
  fEventID = -1;
}

DEFINE_ART_MODULE(extrapolation::PrepNDDeposPackets)

