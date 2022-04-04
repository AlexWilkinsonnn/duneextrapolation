////////////////////////////////////////////////////////////////////////
// Class:       GetNDFDSummedHitMeans
// Plugin Type: analyzer (Unknown Unknown)
// File:        GetNDFDSummedHitMeans_module.cc
//
// Generated Mon Apr 4 2022 by Alexander Wilkinson
//
// Print the mean of the total hit integrals for FD and ND hits.
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
#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <numeric>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TInterpreter.h"

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
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace extrapolation {
  class GetNDFDSummedHitMeans;
}

class extrapolation::GetNDFDSummedHitMeans : public art::EDAnalyzer {
public:
  explicit GetNDFDSummedHitMeans(fhicl::ParameterSet const& p);

  GetNDFDSummedHitMeans(GetNDFDSummedHitMeans const&) = delete;
  GetNDFDSummedHitMeans(GetNDFDSummedHitMeans&&) = delete;
  GetNDFDSummedHitMeans& operator=(GetNDFDSummedHitMeans const&) = delete;
  GetNDFDSummedHitMeans& operator=(GetNDFDSummedHitMeans&&) = delete;

  void analyze(const art::Event& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  std::vector<double> fNDSummedIntegrals;
  std::vector<double> fFDSummedIntegrals;

  std::string fNDPacketLabel;
  std::string fFDHitLabel;
};

extrapolation::GetNDFDSummedHitMeans::GetNDFDSummedHitMeans(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fNDPacketLabel  (p.get<std::string>("NDPacketLabel")),
  fFDHitLabel     (p.get<std::string>("FDHitLabel"))
{
  consumes<std::vector<recob::Hit>>(fNDPacketLabel);
}

void extrapolation::GetNDFDSummedHitMeans::analyze(const art::Event& e)
{
  const auto hitsND = e.getValidHandle<std::vector<recob::Hit>> (fNDPacketLabel);
  const auto hitsFD = e.getValidHandle<std::vector<recob::Hit>> (fFDHitLabel);

  double NDSummedIntegral;
  for (const recob::Hit hit : *hitsND) {
    NDSummedIntegral += hit.Integral();
  }
  fNDSummedIntegrals.push_back(NDSummedIntegral);

  double FDSummedIntegral;
  for (const recob::Hit hit : *hitsFD) {
    FDSummedIntegral += hit.Integral();
  }
  fFDSummedIntegrals.push_back(FDSummedIntegral);
}

void extrapolation::GetNDFDSummedHitMeans::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
  
  std::cout << "Reading ND hits from " << fNDPacketLabel << "\n";
  std::cout << "Reading FD hits from " << fFDHitLabel << "\n";
}

void extrapolation::GetNDFDSummedHitMeans::endJob()
{
  double FDMeanSummedIntegral = std::accumulate(fFDSummedIntegrals.begin(), fFDSummedIntegrals.end(), 0.0); 
  double NDMeanSummedIntegral = std::accumulate(fNDSummedIntegrals.begin(), fNDSummedIntegrals.end(), 0.0); 
  
  std::cout << "#########################################################\n";
  std::cout << "FD mean summed integral = " << FDMeanSummedIntegral << "\n";
  std::cout << "ND mean summed integral = " << NDMeanSummedIntegral << "\n";
  std::cout << "#########################################################\n";
}

DEFINE_ART_MODULE(extrapolation::GetNDFDSummedHitMeans)
