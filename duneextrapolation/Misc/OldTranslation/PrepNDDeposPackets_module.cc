////////////////////////////////////////////////////////////////////////
// Class:       PrepNDDeposPackets
// Plugin Type: producer (Unknown Unknown)
// File:        PrepNDDeposPackets_module.cc
//
// Generated at Thu Apr 07 2022 by Alexander Wilkinson.
//
// Module to read ND packets and deposition from TTree outputtted by
// export_depos_packets_toroot.py. The geometry service is used to
// project packets onto wire + tick and calculate other FD projection
// related quantities, they are written out to a flat tree with an ID
// that is also written into a SED to keep track of. The depositions are
// save as SEDs in the art-root file. All non-active volume (in ND LAr)
// depositiona are projected to a wire + tick and written out along with
// the packets in order to create an infill mask later.
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
#include <set>

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
  std::vector<std::vector<int>>    fInfillMaskZ;
  std::vector<std::vector<int>>    fInfillMaskU;
  std::vector<std::vector<int>>    fInfillMaskV;
  int                              fEventID;
  std::vector<double>              fVertexOut;

  // fhicl params
  std::string  fNDDataLoc;
  unsigned int fCIndex;
  unsigned int fTIndex;
  double       fYShift;
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
  fTreePacketProjections->Branch("infillmaskz", &fInfillMaskZ);
  fTreePacketProjections->Branch("infillmasku", &fInfillMaskU);
  fTreePacketProjections->Branch("infillmaskv", &fInfillMaskV);
  fTreePacketProjections->Branch("eventid", &fEventID);
  fTreePacketProjections->Branch("vertex", &fVertexOut);
}

void extrapolation::PrepNDDeposPackets::produce(art::Event& e)
{
  this->reset();

  if (fEntry < fNEntries) {
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
    // std::cout << detProp.ConvertTicksToX(2000, fPIDZ) << "\n"; = 163.705
    const geo::PlaneGeo pGeoZ = fGeom->Plane(fPIDZ);
    const geo::PlaneGeo pGeoU = fGeom->Plane(fPIDU);
    const geo::PlaneGeo pGeoV = fGeom->Plane(fPIDV);
    const geo::BoxBoundedGeo pGeoBoxZ = pGeoZ.BoundingBox();
    const geo::BoxBoundedGeo pGeoBoxU = pGeoU.BoundingBox();
    const geo::BoxBoundedGeo pGeoBoxV = pGeoV.BoundingBox();

    fTreeDeposPackets->GetEntry(fEntry);

    // Get Vertex alignment info
    double vtxX = fVertex->at(2); // ND z is FD x (drift coordinate)
    double vtxTick = detProp.ConvertXToTicks(vtxX, fPIDZ); // Align wrt 2000 Z plane ticks
    double tickDiff = 2000 - vtxTick;
    double XShift = tickDiff * detProp.GetXTicksCoefficient(fTID.TPC, fTID.Cryostat);

    double vtxZ = fVertex->at(0);
    double ZShift = 0.0;
    if (vtxZ < pGeoBoxZ.MinZ() + 10 || vtxZ > pGeoBoxZ.MaxZ() - 10) { 
      ZShift = (pGeoBoxZ.MinZ() + 50) - vtxZ;
    }

    // std::cout << fGeom->TPC(fTID).MinX() << " - " << fGeom->TPC(fTID).MaxX() << "\n";
    // std::cout << fGeom->TPC(fTID).MinY() << " - " << fGeom->TPC(fTID).MaxY() << "\n";
    // std::cout << fGeom->TPC(fTID).MinZ() << " - " << fGeom->TPC(fTID).MaxZ() << "\n";
    // std::cout << pGeoBoxZ.MinX() << " - " << pGeoBoxZ.MaxX() << ", " << pGeoBoxZ.MinY() << " - " << pGeoBoxZ.MaxY() << ", " << pGeoBoxZ.MinZ() << " - " << pGeoBoxZ.MaxZ() << "\n";
    // std::cout << pGeoBoxU.MinX() << " - " << pGeoBoxU.MaxX() << ", " << pGeoBoxU.MinY() << " - " << pGeoBoxU.MaxY() << ", " << pGeoBoxU.MinZ() << " - " << pGeoBoxU.MaxZ() << "\n";
    // std::cout << pGeoBoxV.MinX() << " - " << pGeoBoxV.MaxX() << ", " << pGeoBoxV.MinY() << " - " << pGeoBoxV.MaxY() << ", " << pGeoBoxV.MinZ() << " - " << pGeoBoxV.MaxZ() << "\n";

    // Make and add the ND depos
    auto SEDs = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
    auto evNum = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

    if (!fOnlyProjections) {
      std::map<geo::View_t, std::set<std::pair<raw::ChannelID_t, unsigned int>>>
        infillMaskChTIcks = { { geo::kZ, { } }, { geo::kU, { } }, { geo::kV, { } } };

      for (const std::vector<double>& depo : *fDepos) {
        int trackID = (int)depo[0];
        int pdg = (int)depo[1];
        double zMin = depo[2] + ZShift;
        double zMax = depo[3] + ZShift;
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

        // A better way to do this is cut out these depos before larnd-sim stage to ensure there
        // no packets partially induced by depos that get cut here.
        try {
          if (fGeom->PositionToTPC(posStart).ID() != fTID) {
            continue;
          }
        }
        catch(...) { // Also want to skip if between TPCs
          continue;
        }

        sim::SimEnergyDeposit SED = sim::SimEnergyDeposit(
            0, electrons, 0, dE, posStart, posEnd, tMin, tMax, trackID, pdg);
        SEDs->push_back(SED);

        double activeVol = depo[12]; // 1.0 if active, -1.0 of not
        if (activeVol < 0) {
          // Didn't need to cut out depos that are goind into SEDs since these will be cut when
          // WC drifts them anyway
          // If deposition is not going to drifted by WC we don't want it
          if (std::abs(posStart.X()) < 5.239 || std::abs(posStart.X()) > 362.916) {
            continue;
          }

          if (pGeoBoxZ.ContainsY(posStart.Y()) &&
              (posStart.Z() > (pGeoBoxZ.MinZ() - (pGeoZ.WirePitch()/2.1)) &&
               posStart.Z() < (pGeoBoxZ.MaxZ() + (pGeoZ.WirePitch()/2.1)))) {
            raw::ChannelID_t ch = fGeom->NearestChannel(posStart, fPIDZ) - fGeom->FirstChannelInROP(fRIDZ);

            double tickRaw = detProp.ConvertXToTicks(xMin, fPIDZ);
            tickRaw += fTickShiftZ;
            unsigned int tick = (unsigned int)tickRaw;

            infillMaskChTIcks[geo::kZ].insert(std::make_pair(ch, tick));
          }

          if (pGeoBoxU.ContainsY(posStart.Y()) &&
              (posStart.Z() > (pGeoBoxU.MinZ() - (pGeoU.WirePitch()/2.1)) &&
               posStart.Z() < (pGeoBoxU.MaxZ() + (pGeoU.WirePitch()/2.1)))) {
            raw::ChannelID_t ch = fGeom->NearestChannel(posStart, fPIDU) - fGeom->FirstChannelInROP(fRIDU);

            double tickRaw = detProp.ConvertXToTicks(xMin, fPIDU);
            tickRaw += fTickShiftU;
            unsigned int tick = (unsigned int)tickRaw;

            infillMaskChTIcks[geo::kU].insert(std::make_pair(ch, tick));
          }

          if (pGeoBoxV.ContainsY(posStart.Y()) &&
              (posStart.Z() > (pGeoBoxV.MinZ() - (pGeoV.WirePitch()/2.1)) &&
               posStart.Z() < (pGeoBoxV.MaxZ() + (pGeoV.WirePitch()/2.1)))) {
            raw::ChannelID_t ch = fGeom->NearestChannel(posStart, fPIDV) - fGeom->FirstChannelInROP(fRIDV);

            double tickRaw = detProp.ConvertXToTicks(xMin, fPIDV);
            tickRaw += fTickShiftV;
            unsigned int tick = (unsigned int)tickRaw;

            infillMaskChTIcks[geo::kV].insert(std::make_pair(ch, tick));
          }
        }
      }

      for (std::pair<raw::ChannelID_t, unsigned int> chTick : infillMaskChTIcks[geo::kZ] ) {
        std::vector<int> infillMask(2, 0);
        infillMask[0] = (int)chTick.first;
        infillMask[1] = (int)chTick.second;
        fInfillMaskZ.push_back(infillMask);
      }

      for (std::pair<raw::ChannelID_t, unsigned int> chTick : infillMaskChTIcks[geo::kU] ) {
        std::vector<int> infillMask(2, 0);
        infillMask[0] = (int)chTick.first;
        infillMask[1] = (int)chTick.second;
        fInfillMaskU.push_back(infillMask);
      }

      for (std::pair<raw::ChannelID_t, unsigned int> chTick : infillMaskChTIcks[geo::kV] ) {
        std::vector<int> infillMask(2, 0);
        infillMask[0] = (int)chTick.first;
        infillMask[1] = (int)chTick.second;
        fInfillMaskV.push_back(infillMask);
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
      double z = packet[0] + ZShift;
      double y = packet[1] + fYShift;
      double x = packet[2] + XShift;
      double adc = packet[4];
      double NDDrift = packet[5];
      double NDModuleX = packet[6]; // x coord (z in FD) relative to ND drift volume

      geo::Point_t packetLoc(x, y, z);

      // If deposition is not going to drifted by WC we don't want it
      if (std::abs(packetLoc.X()) < 5.239 || std::abs(packetLoc.X()) > 362.916) {
        continue;
      }

      try {
        if (fGeom->PositionToTPC(packetLoc).ID() != fTID) {
          continue;
        }
      }
      catch(...) { // Also want to skip if between TPCs
        continue;
      }

      std::vector<double> projection(18, 0.0);
      projection[0] = z;
      projection[1] = y;
      projection[2] = x;
      projection[9] = adc;
      projection[10] = NDDrift;
      projection[17] = NDModuleX;

      // 2.1 to be safe from NearestChannel exceptions
      if (pGeoBoxZ.ContainsY(packetLoc.Y()) &&
          (packetLoc.Z() > (pGeoBoxZ.MinZ() - (pGeoZ.WirePitch()/2.1)) &&
           packetLoc.Z() < (pGeoBoxZ.MaxZ() + (pGeoZ.WirePitch()/2.1)))) {
        raw::ChannelID_t ch = fGeom->NearestChannel(packetLoc, fPIDZ) - fGeom->FirstChannelInROP(fRIDZ);

        double tickRaw = detProp.ConvertXToTicks(x, fPIDZ);
        tickRaw += fTickShiftZ;
        unsigned int tick = (unsigned int)tickRaw;

        double driftDistance = pGeoZ.DistanceFromPlane(packetLoc);

        double wireCoord = pGeoZ.WireCoordinate(packetLoc);
        double wireDistance = (wireCoord - (double)(int)(0.5 + wireCoord)) * pGeoZ.WirePitch();

        projection[3] = ch;
        projection[4] = tick;
        projection[11] = driftDistance;
        projection[14] = wireDistance;
      }
      else {
        projection[3] = -1;
      }

      if (pGeoBoxU.ContainsY(packetLoc.Y()) &&
          (packetLoc.Z() > (pGeoBoxU.MinZ() - (pGeoU.WirePitch()/2.1)) &&
           packetLoc.Z() < (pGeoBoxU.MaxZ() + (pGeoU.WirePitch()/2.1)))) {
        raw::ChannelID_t ch = fGeom->NearestChannel(packetLoc, fPIDU) - fGeom->FirstChannelInROP(fRIDU);

        double tickRaw = detProp.ConvertXToTicks(x, fPIDU);
        tickRaw += fTickShiftU;
        unsigned int tick = (unsigned int)tickRaw;

        double driftDistance = pGeoU.DistanceFromPlane(packetLoc);

        double wireCoord = pGeoU.WireCoordinate(packetLoc);
        double wireDistance = (wireCoord - (double)(int)(0.5 + wireCoord)) * pGeoU.WirePitch();

        projection[5] = ch;
        projection[6] = tick;
        projection[12] = driftDistance;
        projection[15] = wireDistance;
      }
      else {
        projection[5] = -1;
      }

      if (pGeoBoxV.ContainsY(packetLoc.Y()) &&
          (packetLoc.Z() > (pGeoBoxV.MinZ() - (pGeoV.WirePitch()/2.1)) &&
           packetLoc.Z() < (pGeoBoxV.MaxZ() + (pGeoV.WirePitch()/2.1)))) {
        raw::ChannelID_t ch = fGeom->NearestChannel(packetLoc, fPIDV) - fGeom->FirstChannelInROP(fRIDV);

        double tickRaw = detProp.ConvertXToTicks(x, fPIDV);
        tickRaw += fTickShiftV;
        unsigned int tick = (unsigned int)tickRaw;

        double driftDistance = pGeoV.DistanceFromPlane(packetLoc);

        double wireCoord = pGeoV.WireCoordinate(packetLoc);
        double wireDistance = (wireCoord - (double)(int)(0.5 + wireCoord)) * pGeoV.WirePitch();

        projection[7] = ch;
        projection[8] = tick;
        projection[13] = driftDistance;
        projection[16] = wireDistance;
      }
      else {
        projection[7] = -1;
      }

      fPacketProjection.push_back(projection);

      /* // I dont think this will work now that I slightly refactored the way the projection
         // vector gets filled, fix it when I need to.
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
      */
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
  for (geo::PlaneID pID : fGeom->Iterate<geo::PlaneID>(fTID)) {
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

