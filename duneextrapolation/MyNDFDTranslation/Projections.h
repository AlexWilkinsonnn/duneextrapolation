////////////////////////////////////////////////////////////////////////
// Class: Projections
// File:  Projections.h
//
// Created on Tue 21 Jun 22 by Alex Wilkinson.
////////////////////////////////////////////////////////////////////////

#ifndef EXTRAPOLATION_PROJECTIONS_H
#define EXTRAPOLATION_PROJECTIONS_H

#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "duneextrapolation/MyNDFDTranslation/Types.h"

#include <map>
#include <vector>
#include <set>

namespace extrapolation
{

  // Project ND packets onto wire planes and store results.
  class Projections
  {
  public:
    Projections() { fNumPackets = 0; };
    // Constructor that sets projection config
    Projections(const double& tickShiftZ, const double& tickShiftU, const double& tickShiftV,
       const geo::GeometryCore* geom, const bool& calcWireDistance = false);

    // Set the projection config
    void SetProjectionConfig(const double& tickShiftZ, const double& tickShiftU,
       const double& tickShiftV, const geo::GeometryCore* geom,
       const bool& calcWireDistance = false);
    // Passing the new memory address of the new detProp object which varies with event.
    void SetDetProp(const detinfo::DetectorPropertiesData* detProp);

    // Add an ND packet
    void Add(const geo::Point_t& pos, const int& adc, const double& NDDrift);
    void Add(const sim::SimEnergyDeposit& sed);

    // Get number of ND packets stored
    int Size() { return fPackets.size(); }

    // Clear ND packets and projections
    void Clear();

    // Perform the projection
    void ProjectToWires();

    // Get projection information
    std::vector<ProjectionData> GetProjectionData(const readout::ROPID& rID) { return fProjections[rID]; }
    std::vector<readout::ROPID> ActiveROPIDs();
    int GetNumInvalidProjections() { return fNumInvalidProjections; }

  private:
    // Packet data
    std::vector<PacketData> fPackets;

    // Projection config data
    const geo::GeometryCore*               fGeom;
    const detinfo::DetectorPropertiesData* fDetProp;
    // Needed to line up with what WC detsim does
    std::map<geo::View_t, double>          fTickShifts;
    bool                                   fCalcWireDistance;

    // Projection data
    std::map<readout::ROPID, std::vector<ProjectionData>> fProjections;
    int                                                   fNumInvalidProjections;
  };

}

#endif // EXTRAPOLATION_PROJECTIONS_H
