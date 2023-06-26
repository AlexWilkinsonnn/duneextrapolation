////////////////////////////////////////////////////////////////////////
// Class:       GenieTruthDump
// Plugin Type: analyzer (Unknown Unknown)
// File:        GenieTruthDump_module.cc
//
// Crated on 26 June 23 Alex Wilkinson
// Dump genie 4-vector information (position and momentum direction)
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

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <string>
#include <vector>
#include <map>

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"

namespace extrapolation {
  class GenieTruthDump;
}


class extrapolation::GenieTruthDump : public art::EDAnalyzer {
public:
  explicit GenieTruthDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GenieTruthDump(GenieTruthDump const&) = delete;
  GenieTruthDump(GenieTruthDump&&) = delete;
  GenieTruthDump& operator=(GenieTruthDump const&) = delete;
  GenieTruthDump& operator=(GenieTruthDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My function.
  void reset();

private:
  const geo::GeometryCore* fGeom;

  // Labels from fcl
  std::string fMCTruthLabel;

  TTree* fTreeTruth;
  int    fEventNum;
  // MCTruth 4-vector info
  double fX;
  double fY;
  double fZ;
  double fP_x;
  double fP_y;
  double fP_z;
};


extrapolation::GenieTruthDump::GenieTruthDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fMCTruthLabel (p.get<std::string>("MCTruthLabel"))
{
  consumes<std::vector<simb::MCTruth>>(fMCTruthLabel);

  art::ServiceHandle<art::TFileService> tfs;

  fTreeTruth = tfs->make<TTree>("genie", "genie");
  fTreeTruth->Branch("eventnum", &fEventNum, "eventnum/I");
  // MCTruth 4-vector info
  fTreeTruth->Branch("x", &fX, "x/D");
  fTreeTruth->Branch("y", &fY, "y/D");
  fTreeTruth->Branch("z", &fZ, "z/D");
  fTreeTruth->Branch("p_x", &fP_x, "p_x/D");
  fTreeTruth->Branch("p_y", &fP_y, "p_y/D");
  fTreeTruth->Branch("p_z", &fP_z, "p_z/D");
}

void extrapolation::GenieTruthDump::analyze(art::Event const& e)
{
  this->reset();

  fEventNum = e.id().event();

  const auto mcTruth = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);
  const simb::MCParticle nuParticle = mcTruth->at(0).GetNeutrino().Nu();

  fX = nuParticle.Vx();
  fY = nuParticle.Vy();
  fZ = nuParticle.Vz();
  fP_x = nuParticle.Px();
  fP_y = nuParticle.Py();
  fP_z = nuParticle.Pz();

  fTreeTruth->Fill();

  // tan theta = p_x / p_z
  // tan phi = p_y / p_z
  // double p_theta = fP_x > 0.0 ? TMath::ATan(fP_x / fP_z) : -TMath::ATan(fP_x / fP_z);
  // double p_phi = fP_y > 0.0 ? TMath::ATan(fP_y / fP_z) : -TMath::ATan(fP_y / fP_z);
}

void extrapolation::GenieTruthDump::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
}

void extrapolation::GenieTruthDump::endJob()
{
}

void extrapolation::GenieTruthDump::reset()
{
  fEventNum = -99;

  fX = -99.0;
  fY = -99.0;
  fZ = -99.0;
  fP_x = -99.0;
  fP_y = -99.0;
  fP_z = -99.0;
}

DEFINE_ART_MODULE(extrapolation::GenieTruthDump)
