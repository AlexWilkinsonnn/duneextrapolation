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

#include <map>
#include <vector>
#include <set>

extrapolation::Projections::Projections(const double& tickShiftZ, const double& tickShiftU,
  const double& tickShiftV, const geo::GeometryCore* geom,
  const bool& calcWireDistance /* = false */)
{
  SetProjectionConfig(tickShiftZ, tickShiftU, tickShiftV, geom, calcWireDistance);
  fNumPackets = 0;
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
  fPositions.push_back(pos);
  fAdcs.push_back(adc);
  fNDDrifts.push_back(NDDrift);
  fNumPackets++;
}

void extrapolation::Projections::Add(const sim::SimEnergyDeposit& sed)
{
  Add(sed.Start(), sed.NumElectrons(), sed.ScintYieldRatio());
}

int extrapolation::Projections::Size()
{
  return fNumPackets;
}

void extrapolation::Projections::Clear()
{
  fPositions.clear();
  fAdcs.clear();
  fLocalChs.clear();
  fTicks.clear();
  fFDDrifts.clear();
  fWireDistances.clear();
  fActiveROPs.clear();
  fNumPackets = 0;
}

void extrapolation::Projections::ProjectToWires()
{
  // Make vectors of empty maps. Any packets that fall outside drift bb for a plane will have
  // empty map entries.
  fLocalChs = std::vector<std::map<readout::ROPID, int>>(fNumPackets);
  fTicks = std::vector<std::map<readout::ROPID, int>>(fNumPackets);
  fFDDrifts = std::vector<std::map<readout::ROPID, double>>(fNumPackets);
  if (fCalcWireDistance) {
    fWireDistances = std::vector<std::map<readout::ROPID, double>>(fNumPackets);
  }

  for (int i = 0; i < fNumPackets; i++) {
    const geo::Point_t& packetPos = fPositions[i];

    // Check if packet falls outside of WC bb X range. In this case WC wouldnt drift the charge
    // so we wont either for now. Doing this here since it does not depend on plane.
    // TODO Make these number configurable.
    // TODO Confirm these numbers, ~600 too high, the drift is ~[-360, 360].
    // Should only have to worry about the lower bound if xShift is correct since FD is wider than
    // ND.
    if (std::abs(packetPos.X()) < 5.239) { //|| std::abs(packetPos.X()) > 600.019) {
      fNumInvalidProjections += 3;
      continue;
    }

    // Find the 3 wire planes and associated readout planes that the packet will be projected onto
    const geo::TPCID tID = fGeom->PositionToTPCID(packetPos);
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
      if (!pGeoBox.ContainsY(packetPos.Y())) {
        fNumInvalidProjections++;
        continue;
      }
      // This is how WC does it
      if (packetPos.Z() < (pGeoBox.MinZ() - (pGeo.WirePitch()/2))
           || packetPos.Z() > (pGeoBox.MaxZ() + (pGeo.WirePitch()/2))) {
        fNumInvalidProjections++;
        continue;
      }

      // ROP is active if we made it this far
      fActiveROPs.insert(rID);

      // Local channel wrt first channel in ROP ie. column in pixel map
      int ch = (int)(fGeom->NearestChannel(packetPos, pID) - fGeom->FirstChannelInROP(rID));
      fLocalChs[i][rID] = ch;

      // Tick ie. row in pixel map
      std::cout << packetPos.X() << "\n";
      double tickRaw = fDetProp->ConvertXToTicks(packetPos.X(), pID);
      std::cout << tickRaw << "\n";
      tickRaw += fTickShifts[fGeom->View(rID)];
      int tick = (int)tickRaw;
      fTicks[i][rID] = tick;

      // FD drift distance
      double FDDrift = pGeo.DistanceFromPlane(packetPos);
      fFDDrifts[i][rID] = FDDrift;

      if (fCalcWireDistance) {
        double wireCoord = pGeo.WireCoordinate(packetPos);
        double wireDistance = (wireCoord - (double)(int)(0.5 + wireCoord)) * pGeo.WirePitch();
        fWireDistances[i][rID] = wireDistance;
      }
    } // for (int j = 0; j < 3; j++)
  } // for (int i = 0; i < fNumPackets; i++)
}

