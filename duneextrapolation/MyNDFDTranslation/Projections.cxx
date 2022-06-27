////////////////////////////////////////////////////////////////////////
// Class: Projections
// File:  Projections.cxx
//
// Created on Tue 21 Jun 22 by Alex Wilkinson.
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "duneextrapolation/MyNDFDTranslation/Projections.h"
#include "duneextrapolation/MyNDFDTranslation/Types.h"

#include <map>
#include <vector>
#include <set>

extrapolation::Projections::Projections(const double& tickShiftZ, const double& tickShiftU,
  const double& tickShiftV, const geo::GeometryCore* geom,
  const bool& calcWireDistance /* = false */)
{
  SetProjectionConfig(tickShiftZ, tickShiftU, tickShiftV, geom, calcWireDistance);
}

void extrapolation::Projections::SetProjectionConfig(const double& tickShiftZ,
  const double& tickShiftU, const double& tickShiftV, const geo::GeometryCore* geom,
  const bool& calcWireDistance /* = false */)
{
  fTickShifts[geo::kZ] = tickShiftZ;
  fTickShifts[geo::kU] = tickShiftU;
  fTickShifts[geo::kV] = tickShiftV;
  fGeom = geom;
  fCalcWireDistance = calcWireDistance;
}

void extrapolation::Projections::SetDetProp(const detinfo::DetectorPropertiesData* detProp)
{
  fDetProp = detProp;
}

void extrapolation::Projections::Add(const geo::Point_t& pos, const int& adc, const double& NDDrift)
{
  fPackets.push_back(PacketData(pos, adc, NDDrift));
}

void extrapolation::Projections::Add(const sim::SimEnergyDeposit& sed)
{
  Add(sed.Start(), sed.NumElectrons(), sed.ScintYieldRatio());
}

void extrapolation::Projections::Clear()
{
  fPackets.clear();
  fProjections.clear();
}

void extrapolation::Projections::ProjectToWires()
{
  for (PacketData& packet : fPackets) {
    // Check if packet falls outside of WC bb X range. In this case WC wouldnt drift the charge
    // so we wont either for now. Doing this here since it does not depend on plane.
    // TODO Make these number configurable.
    // TODO Confirm these numbers, ~600 too high, the drift is ~[-360, 360].
    // Should only have to worry about the lower bound if xShift is correct since FD is wider than
    // ND.
    if (std::abs(packet.pos.X()) < 5.239) { //|| std::abs(packet.pos.X()) > 600.019) {
      fNumInvalidProjections += 3;
      continue;
    }

    // Find the 3 wire planes and associated readout planes that the packet will be projected onto
    const geo::TPCID tID = fGeom->PositionToTPCID(packet.pos);
    std::vector<geo::PlaneID> pIDs;
    std::vector<readout::ROPID> rIDs;

    for(geo::PlaneID pID : fGeom->IteratePlaneIDs(tID)){
      if (std::count(pIDs.begin(), pIDs.end(), pID) || std::count(rIDs.begin(), rIDs.end(), fGeom->WirePlaneToROP(pID))) {
        std::cout << "!!!!\nconfusion over iterateplaneids in ProjectToWires\n!!!!\n";
      }
      pIDs.push_back(pID);
      rIDs.push_back(fGeom->WirePlaneToROP(pID));
    }
    if (pIDs.size() != 3) {
      std::cout << "!!!!\nconfusion over iterateplaneids in ProjectToWires\n!!!!\n";
    }

    // Project to each of the 3 wire planes
    for (int j = 0; j < 3; j++) {
      const geo::PlaneID& pID = pIDs[j];
      const geo::PlaneGeo pGeo = fGeom->Plane(pID);
      const readout::ROPID& rID = rIDs[j];

      // Check if packet has valid Y Z to be drifted
      const geo::BoxBoundedGeo pGeoBox = pGeo.BoundingBox();

      // This should never happen since ND shorter than FD but good to check
      if (!pGeoBox.ContainsY(packet.pos.Y())) {
        fNumInvalidProjections++;
        continue;
      }
      // This is how WC does it
      if (packet.pos.Z() < (pGeoBox.MinZ() - (pGeo.WirePitch()/2))
           || packet.pos.Z() > (pGeoBox.MaxZ() + (pGeo.WirePitch()/2))) {
        fNumInvalidProjections++;
        continue;
      }

      ProjectionData projData(packet);

      // Local channel wrt first channel in ROP ie. column in pixel map
      int ch = (int)(fGeom->NearestChannel(packet.pos, pID) - fGeom->FirstChannelInROP(rID));
      projData.localCh = ch;

      // Tick ie. row in pixel map
      double tickRaw = fDetProp->ConvertXToTicks(packet.pos.X(), pID);
      tickRaw += fTickShifts[fGeom->View(rID)];
      int tick = (int)tickRaw;
      projData.tick = tick;

      // FD drift distance
      double FDDrift = pGeo.DistanceFromPlane(packet.pos);
      projData.FDDrift = FDDrift;

      if (fCalcWireDistance) {
        double wireCoord = pGeo.WireCoordinate(packet.pos);
        double wireDistance = (wireCoord - (double)(int)(0.5 + wireCoord)) * pGeo.WirePitch();
        projData.wireDistance = wireDistance;
      }

      fProjections[rID].push_back(projData);
      fActiveChTicks[rID].insert(std::make_pair(ch, tick));
    } // for (int j = 0; j < 3; j++)
  } // for (PacketData& packet : fPackets)
}

std::vector<readout::ROPID> extrapolation::Projections::ActiveROPIDs()
{
  std::vector<readout::ROPID> rIDs;

  for (std::pair<readout::ROPID, std::vector<ProjectionData>> proj : fProjections) {
    rIDs.push_back(proj.first);
  }

  return rIDs;
}

