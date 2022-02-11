////////////////////////////////////////////////////////////////////////
// Class:       AddNDDepos
// Plugin Type: producer (Unknown Unknown)
// File:        AddNDDepos_module.cc
//
// Generated at Wed Jan 19 08:14:24 2022 by Alexander Wilkinson using cetskelgen
// from  version .
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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace extrapolation {
  class AddNDDepos;
}

class extrapolation::AddNDDepos : public art::EDProducer {
public:
  explicit AddNDDepos(fhicl::ParameterSet const& p);

  AddNDDepos(AddNDDepos const&) = delete;
  AddNDDepos(AddNDDepos&&) = delete;
  AddNDDepos& operator=(AddNDDepos const&) = delete;
  AddNDDepos& operator=(AddNDDepos&&) = delete;

  void produce(art::Event& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  int fEventNumber;
  int fFileCnt;
  std::string fNDDeposLoc;
  double fYShift;
  double fZShift;
  geo::TPCID fTID;
  geo::PlaneID fPID;
};

extrapolation::AddNDDepos::AddNDDepos(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fNDDeposLoc  (p.get<std::string>("NDDeposLoc"))
{
  produces<std::vector<sim::SimEnergyDeposit>>();
  produces<std::vector<sim::SimEnergyDeposit>>("EventNumber");
}

void extrapolation::AddNDDepos::produce(art::Event& e)
{
  // auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  bool fileExists = false;
  while (!fileExists) {
    std::ifstream iFile;
    std::string fileName = fNDDeposLoc + "/ND_depos_" + std::to_string(fEventNumber) + ".txt";
    iFile.open(fileName);

    if (iFile) {
      fileExists = true;

      // Make and add the sim depos
      auto SEDs = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
      auto evNum = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

      std::string line;
      std::getline(iFile, line); // metadata line

      std::string lineSplitLeft = line.substr(line.find("vtx_z:") + 6);
      std::string lineSplitLeftRight = lineSplitLeft.substr(0, lineSplitLeft.find(','));
      double vtxX = std::stod(lineSplitLeftRight);

      double vtxTick = detProp.ConvertXToTicks(vtxX, fPID) - 7.8; // empirical wirecell offset
      // std::cout << "vtxX=" << vtxX << " vtxTick=" << vtxTick;
      double tickDiff = 2000 - vtxTick;
      double XShift = tickDiff * detProp.GetXTicksCoefficient(fTID.TPC, fTID.Cryostat);
      // std::cout << " XShift=" << XShift << " detProp.ConvertXToTicks(vtxX + XShift, fPID) - 7.8=" << detProp.ConvertXToTicks(vtxX + XShift, fPID) - 7.8 << "\n";

      while (std::getline(iFile, line)) {
        std::istringstream iss(line);

        // X <--> Z for ND and FD. Drift coordinate is X for FD and Z for ND.
        std::string entry;
        std::getline(iss, entry, ',');
        int trackID = std::stoi(entry);
        std::getline(iss, entry, ',');
        int pdg = std::stoi(entry);
        std::getline(iss, entry, ',');
        double zMin = std::stod(entry) + fZShift;
        std::getline(iss, entry, ',');
        double zMax = std::stod(entry) + fZShift;
        std::getline(iss, entry, ',');
        double yMin = std::stod(entry) + fYShift;
        std::getline(iss, entry, ',');
        double yMax = std::stod(entry) + fYShift;
        std::getline(iss, entry, ',');
        double xMin = std::stod(entry) + XShift;
        std::getline(iss, entry, ',');
        double xMax = std::stod(entry) + XShift;
        std::getline(iss, entry, ',');
        double tMin = std::stod(entry);
        std::getline(iss, entry, ',');
        double tMax = std::stod(entry);
        std::getline(iss, entry, ',');
        int electrons = std::stoi(entry);
        std::getline(iss, entry, ',');
        double dE = std::stod(entry);

        // std::cout << trackID << " " << pdg << " ([" << xMin << "," << xMax << "],[" << yMin << "," <<
        //   yMax << "],[" << zMin << "," << zMax << "],[" << tMin << "," << tMax << "]) " << electrons <<
        //   " " << dE << "\n";

        geo::Point_t posStart = geo::Point_t(xMin, yMin, zMin);
        geo::Point_t posEnd = geo::Point_t(xMax, yMax, zMax);

        sim::SimEnergyDeposit SED = sim::SimEnergyDeposit(0, electrons, 0, dE, posStart, posEnd,
          tMin, tMax, trackID, pdg);
        SEDs->push_back(SED);
      }

      // Check
      for (const sim::SimEnergyDeposit& SED : *SEDs) {
        const geo::TPCID tpcID = fGeom->PositionToTPCID(SED.MidPoint());
        if (tpcID != fTID) {
          std::cout << tpcID.TPC << "\n";
        }
      }

      geo::Point_t posStart = geo::Point_t(0, 0, 0);
      geo::Point_t posEnd = geo::Point_t(0, 0, 0);
      sim::SimEnergyDeposit ID = sim::SimEnergyDeposit(0, 0, 0, 0, posStart, posEnd, 0, 0, fEventNumber);
      evNum->push_back(ID);
      fEventNumber++;

      e.put(std::move(SEDs));
      e.put(std::move(evNum), "EventNumber");
    }
    else {
      std::cout << fEventNumber << " " << fFileCnt << "\n";
      if (fEventNumber >= fFileCnt) {
        std::cout << "############\nNumber of events being processed must be less than " << fFileCnt << 
          "\n############\n";
        std::cout << "Too many events for files in input directory!" << std::endl;
      }

      fEventNumber++;
    }
  }
}

void extrapolation::AddNDDepos::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  // Check APAs to choose one
  // for (const geo::PlaneGeo& pGeo : fGeom->IteratePlanes()) {
  //   if (fGeom->View(pGeo.ID()) == geo::kZ) {
  //     std::cout << pGeo.ID() << "\n";
  //     std::cout << pGeo.FirstWire().GetStart().Z() << " " << pGeo.LastWire().GetStart().Z() << "\n";

  //     const geo::BoxBoundedGeo pGeoBox = pGeo.BoundingBox();
  //     std::cout << "z: [" << pGeoBox.MinZ() - (pGeo.WirePitch()/2) << "," <<
  //       (pGeoBox.MaxZ() + (pGeo.WirePitch()/2)) << "]\n";
  //     std::cout << "x: [" << pGeoBox.MinX()<< "," << pGeoBox.MaxX() << "]\n";
  //     std::cout << "y: [" << pGeoBox.MinY()<< "," << pGeoBox.MaxY() << "]\n";
  //   }
  // }

  std::filesystem::path path = fNDDeposLoc;
  fFileCnt = 0;
  for (auto& dirEntry : std::filesystem::directory_iterator(path)) {
    (void)dirEntry;
    fFileCnt++;
  }
  std::cout << "############\nNumber of events being processed must be less than " << fFileCnt << 
    "\n############\n";

  // Chosen wire plane and shifts required to bring Y,Z in front of it
  // X shift will be done using the vertex
  const geo::CryostatID cID(0);
  fTID = geo::TPCID(cID, 10);
  fPID = geo::PlaneID(fTID, 2);
  fYShift = 200; // cm
  fZShift = -14.8615; // cm 

  fEventNumber = 0;
}

void extrapolation::AddNDDepos::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::AddNDDepos)
