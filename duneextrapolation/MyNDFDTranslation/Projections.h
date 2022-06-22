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
    int Size();

    // Clear ND packets and projections
    void Clear();

    // Perform the projection
    void ProjectToWires();

    // Get projection information
    std::set<readout::ROPID> ActiveROPIDs() { return fActiveROPs; }
    std::map<readout::ROPID, int> GetChs(const int& i) { return fLocalChs[i]; }
    std::map<readout::ROPID, int> GetTicks(const int& i) { return fTicks[i]; }
    int GetAdc(const int& i) { return fAdcs[i]; }
    double  GetNDDrift(const int& i) { return fNDDrifts[i]; }
    std::map<readout::ROPID, double> GetFDDrifts(const int& i) { return fFDDrifts[i]; }
    std::map<readout::ROPID, double> GetWireDistances(const int& i) { return fWireDistances[i]; }
    int GetNumInvalidProjections() { return fNumInvalidProjections; }

  private:
    // Packet data
    std::vector<geo::Point_t> fPositions;
    std::vector<int>          fAdcs;
    std::vector<double>       fNDDrifts;
    int                       fNumPackets;

    // Projection config data
    const geo::GeometryCore*               fGeom;
    const detinfo::DetectorPropertiesData* fDetProp;
    // Needed to line up with what WC detsim does
    std::map<geo::View_t, double>          fTickShifts;
    bool                                   fCalcWireDistance;

    // Projection data
    std::vector<std::map<readout::ROPID, int>>    fLocalChs;
    std::vector<std::map<readout::ROPID, int>>    fTicks;
    std::vector<std::map<readout::ROPID, double>> fFDDrifts;
    std::vector<std::map<readout::ROPID, double>> fWireDistances;
    std::set<readout::ROPID>                      fActiveROPs;
    int                                           fNumInvalidProjections;
  };
}

#endif // EXTRAPOLATION_PROJECTIONS_H
