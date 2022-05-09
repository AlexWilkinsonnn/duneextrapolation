////////////////////////////////////////////////////////////////////////
// Class:       Interactive
// Plugin Type: analyzer (Unknown Unknown)
// File:        Interactive_module.cc
//
// Generated Thu Mar 24 2022 by Alexander Wilkinson
//
// A module that just prints crap I care about out
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

namespace interactive {
  class Session;
}

class interactive::Session : public art::EDAnalyzer {
public:
  explicit Session(fhicl::ParameterSet const& p);

  Session(Session const&) = delete;
  Session(Session&&) = delete;
  Session& operator=(Session const&) = delete;
  Session& operator=(Session&&) = delete;

  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  std::string fSEDLabel;
};

interactive::Session::Session(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fSEDLabel (p.get<std::string>("SEDLabel"))
{
  if (fSEDLabel != "none") {
    consumes<std::vector<sim::SimEnergyDeposit>>(fSEDLabel);
  }
}

void interactive::Session::analyze(art::Event const& e)
{
  if (fSEDLabel != "none") {
    const auto SEDs = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fSEDLabel);

    std::vector<geo::Length_t> stepSizes;
    std::vector<int> numElectrons;
    std::vector<double> energies;
    for (auto& SED : *SEDs) {
      stepSizes.push_back(SED.StepLength());
      numElectrons.push_back(SED.NumElectrons());
      energies.push_back(SED.Energy());
    }

    { using std::cout;
      cout << "StepSizes:\n{ ";
      for (auto stepSize : stepSizes) {
        cout << stepSize << ", ";
      }
      cout << " }\n";
      cout << "Min = " << *std::min_element(stepSizes.begin(), stepSizes.end());
      cout << "  Max = " << *std::max_element(stepSizes.begin(), stepSizes.end());
      cout << "  Mean = " << std::accumulate(stepSizes.begin(), stepSizes.end(), 0.0)/stepSizes.size();
      cout << "\n\n";

      cout << "NumElectrons:\n{ ";
      for (auto numElectron : numElectrons) {
        cout << numElectron << ", ";
      }
      cout << " }\n";
      cout << "Min = " << *std::min_element(numElectrons.begin(), numElectrons.end());
      cout << "  Max = " << *std::max_element(numElectrons.begin(), numElectrons.end());
      cout << "  Mean = " << std::accumulate(numElectrons.begin(), numElectrons.end(), 0.0)/numElectrons.size();
      cout << "\n\n";

      // cout << "Energies:\n{ ";
      // for (auto energy : energies) { //   cout << energy << ", ";
      // }
      // cout << " }\n";
      // cout << "Min = " << *std::min_element(energies.begin(), energies.end());
      // cout << "  Max = " << *std::max_element(energies.begin(), energies.end());
      // cout << "  Mean = " << std::accumulate(energies.begin(), energies.end(), 0.0)/energies.size();
      // cout << "\n\n";
    }
  }
}

void interactive::Session::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  // Examine readout plane I have chosen for nd->fd translation
  if (true) {
    using std::cout;

    unsigned int cryoIndex = 0;
    unsigned int tpcIndex = 10;
    unsigned int planeIndexZ = 2;
    unsigned int planeIndexU = 1;
    unsigned int planeIndexV = 0;

    const geo::CryostatID cID(cryoIndex);
    const geo::TPCID tID(cID, tpcIndex);
    const geo::PlaneID pIDZ(tID, planeIndexZ);
    const readout::ROPID rIDZ = fGeom->WirePlaneToROP(pIDZ);
    const geo::PlaneID pIDU(tID, planeIndexU);
    const readout::ROPID rIDU = fGeom->WirePlaneToROP(pIDU);
    const geo::PlaneID pIDV(tID, planeIndexV);
    const readout::ROPID rIDV = fGeom->WirePlaneToROP(pIDV);

    cout << "FirstChannelInROP Z = " << fGeom->FirstChannelInROP(rIDZ) << "\n";
    if (fGeom->View(rIDU) == geo::kU) { // Not sure if I have plane ID numbers the right way round
      cout << "FirstChannelInROP U = " << fGeom->FirstChannelInROP(rIDU) << "\n";
      cout << "FirstChannelInROP V = " << fGeom->FirstChannelInROP(rIDV) << "\n";
    }
    else {
      cout << "FirstChannelInROP U = " << fGeom->FirstChannelInROP(rIDV) << "\n";
      cout << "FirstChannelInROP V = " << fGeom->FirstChannelInROP(rIDU) << "\n";
    }

    const geo::TPCGeo tGeo = fGeom->TPC(tID);
    const geo::BoxBoundedGeo tBBGeo = tGeo.BoundingBox();
    cout << "tBBGeo.MinX()=" << tBBGeo.MinX() << ", tBBGeo.MaxX()=" << tBBGeo.MaxX() << "\n";
  } 

  // Examine induction plane I will use
  if (false) {
    using std::cout;

    unsigned int cryoIndex = 0;
    unsigned int tpcIndex = 10;
    unsigned int planeIndex = 1;

    const geo::CryostatID cID(cryoIndex);
    const geo::TPCID tID(cID, tpcIndex);
    const geo::PlaneID pID(tID, planeIndex);

    const geo::PlaneGeo pGeo = fGeom->Plane(pID);

    cout << "ThetaZ = " << pGeo.ThetaZ() << "\n";
    cout << "FirstWire Start Y = " << pGeo.FirstWire().GetStart().Y() << "\n";
    cout << "LastWire Start Y = " << pGeo.LastWire().GetStart().Y() << "\n";
    cout << "LastWire ThetaZ = " << pGeo.LastWire().ThetaZ() << "\n";
  }

  if (false) {
    using std::cout;

    std::cout << "Z wire pitch = " << fGeom->WirePitch(geo::kZ) <<
      ", U wire pitch = " << fGeom->WirePitch(geo::kU) <<
      ", V wire pitch = " << fGeom->WirePitch(geo::kV) << "\n";
  }
}

void interactive::Session::endJob()
{
}

DEFINE_ART_MODULE(interactive::Session)
