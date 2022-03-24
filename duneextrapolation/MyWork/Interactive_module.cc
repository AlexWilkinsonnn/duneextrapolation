////////////////////////////////////////////////////////////////////////
// Class:       Interactive
// Plugin Type: analyzer (Unknown Unknown)
// File:        Interactive_module.cc
//
// Generated Thu Mar 24 2022 by Alexander Wilkinson
//
// A module that just uses beginJob() to print crap I care about out
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

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

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
};


interactive::Session::Session(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{

}

void interactive::Session::analyze(art::Event const& e)
{

}

void interactive::Session::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  // Examine readout plane I have chosen for nd->fd translation
  unsigned int CryoIndex = 0;
  unsigned int TpcIndex = 10;
  unsigned int PlaneIndex = 2;

  const geo::CryostatID cID(CryoIndex);
  const geo::TPCID tID(cID, TpcIndex);
  const geo::PlaneID pID(tID, PlaneIndex);
  const readout::ROPID rID = fGeom->WirePlaneToROP(pID);

  std::cout << "FirstChannelInROP = " << fGeom->FirstChannelInROP(rID) << "\n";
}

void interactive::Session::endJob()
{
}

DEFINE_ART_MODULE(interactive::Session)
