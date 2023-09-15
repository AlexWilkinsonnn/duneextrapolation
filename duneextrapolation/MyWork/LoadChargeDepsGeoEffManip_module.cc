////////////////////////////////////////////////////////////////////////
// Class:       LoadChargeDepsGeoEffManip
// Plugin Type: producer (Unknown Unknown)
// File:        LoadChargeDepsGeoEffManip_module.cc
//
// Generated at Sep 13 2023 by Alexander Wilkinson.
//
// Read charge deposition information from another source into a
// SimEnergyDeposit product. Expects the fhicl to point to a root with
// a TTree called myEvents that has the branches:
// float nEdeps[eventSize]
// float deps_E_MeV[eventSize]
// float deps_start_t_us[eventSize]
// float deps_stop_t_us[eventSize]
// float throwVtx_fd cm[3]
// float fdthrow_deps_start_x_cm[eventSize]
// float fdthrow_deps_stop_x_cm[eventSize]
// float fdthrow_deps_start_y_cm[eventSize]
// float fdthrow_deps_stop_y_cm[eventSize]
// float fdthrow_deps_start_z_cm[eventSize]
// float fdthrow_deps_stop_z_cm[eventSize]

// IMPORTANT Above coordinates are in the FD convention where X is the drift
// direction and Z is beam direction. ND convention is to swap these.
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

#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace extrapolation {
  class LoadChargeDepsGeoEffManip;
}

class extrapolation::LoadChargeDepsGeoEffManip : public art::EDProducer {
public:
  explicit LoadChargeDepsGeoEffManip(fhicl::ParameterSet const& p);

  LoadChargeDepsGeoEffManip(LoadChargeDepsGeoEffManip const&) = delete;
  LoadChargeDepsGeoEffManip(LoadChargeDepsGeoEffManip&&) = delete;
  LoadChargeDepsGeoEffManip& operator=(LoadChargeDepsGeoEffManip const&) = delete;
  LoadChargeDepsGeoEffManip& operator=(LoadChargeDepsGeoEffManip&&) = delete;

  void produce(art::Event& e) override;

  void beginJob() override;
  void endJob() override;

  void reset();

private:
  const geo::GeometryCore* fGeom;

  // Reading input tree
  int    fNEntries;
  int    fEntry;
  TTree* fTreeDepos;

  int   fNDeps;
  float fDepEs[100000];
  float fDepStartTs[100000];
  float fDepStopTs[100000];
  float fDepStartXs[100000];
  float fDepStopXs[100000];
  float fDepStartYs[100000];
  float fDepStopYs[100000];
  float fDepStartZs[100000];
  float fDepStopZs[100000];

  // fhicl params
  std::string fDepoDataLoc;
};

extrapolation::LoadChargeDepsGeoEffManip::LoadChargeDepsGeoEffManip(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fDepoDataLoc (p.get<std::string>("DepoDataLoc"))
{
  produces<std::vector<sim::SimEnergyDeposit>>();
}

void extrapolation::LoadChargeDepsGeoEffManip::produce(art::Event& e)
{
  if (fEntry >= fNEntries) {
    std::cout << "Gone beyond number of entries in tree (" << fNEntries << ")\n";
    return;
  }

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  fTreeDepos->GetEntry(fEntry);

  if (fNDeps > 100000) {
    std::cout << "!!!! Number of depositions is too large to fit in array ("
              << fNDeps << "). Prepare for a segfault "
              << "(don't want to skip this event as it will make nd-fd pair matching "
              << "hard currently !!!!\n"; 
  }

  // Make and add the ND depos
  auto SEDs = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

  for (int iDep = 0; iDep < fNDeps; iDep++) { 
    geo::Point_t posStart(fDepStartXs[iDep], fDepStartYs[iDep], fDepStartZs[iDep]);
    geo::Point_t posEnd(fDepStopXs[iDep], fDepStopYs[iDep], fDepStopZs[iDep]);
    sim::SimEnergyDeposit SED = sim::SimEnergyDeposit(
      0, 0, 0, fDepEs[iDep], posStart, posEnd, fDepStartTs[iDep], fDepStopTs[iDep], 1, 1
    );
    SEDs->push_back(SED);
  }

  e.put(std::move(SEDs));

  fEntry++;
}

void extrapolation::LoadChargeDepsGeoEffManip::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  std::cout << "Reading file from " << fDepoDataLoc << "\n";
  TFile* fileDepos = new TFile(fDepoDataLoc.c_str());
  fTreeDepos = (TTree*)fileDepos->Get("myEvents");
  fTreeDepos->SetBranchAddress("nEdeps", &fNDeps);
  fTreeDepos->SetBranchAddress("deps_E_MeV", fDepEs);
  fTreeDepos->SetBranchAddress("fdthrow_deps_start_x_cm", fDepStartXs);
  fTreeDepos->SetBranchAddress("fdthrow_deps_stop_x_cm", fDepStopXs);
  fTreeDepos->SetBranchAddress("fdthrow_deps_start_y_cm", fDepStartYs);
  fTreeDepos->SetBranchAddress("fdthrow_deps_stop_y_cm", fDepStopYs);
  fTreeDepos->SetBranchAddress("fdthrow_deps_start_z_cm", fDepStartZs);
  fTreeDepos->SetBranchAddress("fdthrow_deps_stop_z_cm", fDepStopZs);

  fEntry = 0;
  fNEntries = fTreeDepos->GetEntries();
  std::cout << "File has " << fNEntries << " entries\n";
}


void extrapolation::LoadChargeDepsGeoEffManip::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::LoadChargeDepsGeoEffManip)

