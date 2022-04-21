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
  consumes<std::vector<sim::SimEnergyDeposit>>(fSEDLabel);
}

void interactive::Session::analyze(art::Event const& e)
{
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
    // for (auto energy : energies) {
    //   cout << energy << ", ";
    // }
    // cout << " }\n";
    // cout << "Min = " << *std::min_element(energies.begin(), energies.end());
    // cout << "  Max = " << *std::max_element(energies.begin(), energies.end());
    // cout << "  Mean = " << std::accumulate(energies.begin(), energies.end(), 0.0)/energies.size();
    // cout << "\n\n";
  }
}

void interactive::Session::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  // Examine readout plane I have chosen for nd->fd translation
  if (false) {
    unsigned int CryoIndex = 0;
    unsigned int TpcIndex = 10;
    unsigned int PlaneIndex = 2;

    const geo::CryostatID cID(CryoIndex);
    const geo::TPCID tID(cID, TpcIndex);
    const geo::PlaneID pID(tID, PlaneIndex);
    const readout::ROPID rID = fGeom->WirePlaneToROP(pID);

    std::cout << "FirstChannelInROP = " << fGeom->FirstChannelInROP(rID) << "\n";

    const geo::TPCGeo tGeo = fGeom->TPC(tID);
    const geo::BoxBoundedGeo tBBGeo = tGeo.BoundingBox();
    std::cout << "tBBGeo.MinX()=" << tBBGeo.MinX() << ", tBBGeo.MaxX()=" << tBBGeo.MaxX() << "\n";
  } 

  // Examine induction plane I will use
  if (false) {
    using std::cout;

    unsigned int CryoIndex = 0;
    unsigned int TpcIndex = 10;
    unsigned int PlaneIndex = 1;

    const geo::CryostatID cID(CryoIndex);
    const geo::TPCID tID(cID, TpcIndex);
    const geo::PlaneID pID(tID, PlaneIndex);

    const geo::PlaneGeo pGeo = fGeom->Plane(pID);

    cout << "ThetaZ = " << pGeo.ThetaZ() << "\n";
    cout << "FirstWire Start Y = " << pGeo.FirstWire().GetStart().Y() << "\n";
    cout << "LastWire Start Y = " << pGeo.LastWire().GetStart().Y() << "\n";
    cout << "LastWire ThetaZ = " << pGeo.LastWire().ThetaZ() << "\n";
  }
}

void interactive::Session::endJob()
{
}

DEFINE_ART_MODULE(interactive::Session)
