////////////////////////////////////////////////////////////////////////
// Class:       CheckSpSolveDisambig
// Plugin Type: analyzer (Unknown Unknown)
// File:        CheckSpSolveDisambig_module.cc
//
// Generated at Thu Jan 20 08:15:21 2022 by Alexander Wilkinson using cetskelgen
// from  version .
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
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <vector>
#include <set>
#include <map>


namespace extrapolation {
  class CheckSpSolveDisambig;
}


class extrapolation::CheckSpSolveDisambig : public art::EDAnalyzer {
public:
  explicit CheckSpSolveDisambig(fhicl::ParameterSet const& p);

  CheckSpSolveDisambig(CheckSpSolveDisambig const&) = delete;
  CheckSpSolveDisambig(CheckSpSolveDisambig&&) = delete;
  CheckSpSolveDisambig& operator=(CheckSpSolveDisambig const&) = delete;
  CheckSpSolveDisambig& operator=(CheckSpSolveDisambig&&) = delete;

  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

  void reset();
private:
  const geo::GeometryCore* fGeom;

  std::string fBasicDisambigHitLabel;
  std::string fSpSolveDisambigHitLabel;
};


extrapolation::CheckSpSolveDisambig::CheckSpSolveDisambig(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fBasicDisambigHitLabel   (p.get<std::string>("BasicDisambigHitLabel")),
    fSpSolveDisambigHitLabel (p.get<std::string>("SpSolveDisambigHitLabel"))
{
  consumes<std::vector<recob::Hit>>(fBasicDisambigHitLabel);
  consumes<std::vector<recob::Hit>>(fSpSolveDisambigHitLabel);
}


void extrapolation::CheckSpSolveDisambig::analyze(art::Event const& e)
{
  art::Handle<std::vector<recob::Hit>> basicHits;
  e.getByLabel(fBasicDisambigHitLabel, basicHits);
  art::Handle<std::vector<recob::Hit>> spsolveHits;
  e.getByLabel(fSpSolveDisambigHitLabel, spsolveHits);

  std::cout << "basic: " << basicHits->size() << " hits -- "
            << "spsolve: " << spsolveHits->size() << "hits\n";

  float basicSummedIntegral = 0.0;
  for (recob::Hit hit : *basicHits)
    basicSummedIntegral += hit.Integral();
  float spsolveSummedIntegral = 0.0;
  for (recob::Hit hit : *spsolveHits)
    spsolveSummedIntegral += hit.Integral();
  std::cout << "SummedIntegral:\n" <<  "basic: " << basicSummedIntegral << " -- "
            << "spsolve: " << spsolveSummedIntegral << "\n";

  std::map<unsigned int, float> basicSummedIntegralByWire;
  std::map<unsigned int, std::vector<recob::Hit>> basicHitsByWire;
  for (recob::Hit hit : *basicHits) {
    unsigned int globalIndex = 
      hit.WireID().TPC * 100000 + hit.WireID().Plane * 10000 + hit.WireID().Wire;
    basicSummedIntegralByWire[globalIndex] += hit.Integral();
    basicHitsByWire[globalIndex].push_back(hit);
  }
  std::map<unsigned int, float> spsolveSummedIntegralByWire;
  std::map<unsigned int, std::vector<recob::Hit>> spsolveHitsByWire;
  for (recob::Hit hit : *spsolveHits) { 
    unsigned int globalIndex = 
      hit.WireID().TPC * 100000 + hit.WireID().Plane * 10000 + hit.WireID().Wire;
    spsolveSummedIntegralByWire[globalIndex] += hit.Integral();
    spsolveHitsByWire[globalIndex].push_back(hit);
  }
  std::set<unsigned int> allWires;
  for (std::pair<unsigned int, float> pair : basicSummedIntegralByWire)
    allWires.insert(pair.first);
  for (std::pair<unsigned int, float> pair : spsolveSummedIntegralByWire)
    allWires.insert(pair.first);
  std::cout << "Non-matching SummedIntegral by wire (basic -- spsolve):\n";
  for (unsigned int wire : allWires) { 
    if ((int)basicSummedIntegralByWire[wire] == (int)spsolveSummedIntegralByWire[wire])
      continue;
    std::cout << "T*100000 + P*10000 + W = " << wire << ": "
              << basicSummedIntegralByWire[wire] << " -- "
              << spsolveSummedIntegralByWire[wire] << "\n";
    if (basicSummedIntegralByWire[wire] == 0.0) {
      for (recob::Hit hit : spsolveHitsByWire[wire]) {
        std::cout << hit.StartTick() << ", " << hit.EndTick() << ", "
                  << hit.PeakAmplitude() << ", " << "\n";
      }
    }
  }
}


void extrapolation::CheckSpSolveDisambig::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
}


void extrapolation::CheckSpSolveDisambig::endJob()
{
}


DEFINE_ART_MODULE(extrapolation::CheckSpSolveDisambig)

