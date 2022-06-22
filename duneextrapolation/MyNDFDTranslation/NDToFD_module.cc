////////////////////////////////////////////////////////////////////////
// Class:       NDToFD
// Plugin Type: producer (Unknown Unknown)
// File:        NDToFD_module.cc
//
// Created on Tue 21 Jun 22 by Alex Wilkinson.
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

#include <torch/torch.h>
#include <torch/script.h>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "duneextrapolation/MyNDFDTranslation/Projections.h"

namespace extrapolation {
  class NDToFD;
}

class extrapolation::NDToFD : public art::EDProducer {
public:
  explicit NDToFD(fhicl::ParameterSet const& p);

  NDToFD(NDToFD const&) = delete;
  NDToFD(NDToFD&&) = delete;
  NDToFD& operator=(NDToFD const&) = delete;
  NDToFD& operator=(NDToFD&&) = delete;

  void produce(art::Event& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  std::map<geo::View_t, torch::jit::script::Module> fNetworks;

  Projections fProj;

  // fhicl params
  std::string fNetworkPathZ;
  std::string fNetworkPathU;
  std::string fNetworkPathV;
  int         fPixelMapMode;
  // TODO add an option to save the pixel map to the event. Will need to define my own data product
  // in order to do this.
};

extrapolation::NDToFD::NDToFD(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fNetworkPathZ (p.get<std::string>("NetworkPathZ")),
    fNetworkPathU (p.get<std::string>("NetworkPathU")),
    fNetworkPathV (p.get<std::string>("NetworkPathV")),
    fPixelMapMode (p.get<int>("PixelMapMode"))
{
  // produces<std::vector<raw::RawDigit>>("NDTranslated");
}

void extrapolation::NDToFD::produce(art::Event& e)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
  fProj.SetDetProp(&detProp);

  geo::PlaneID pID(0, 0, 1);
  std::cout << "hello\n";
  std::cout << detProp.Efield() << "\n";
  std::cout << detProp.ConvertXToTicks(100, pID) << "\n";
  std::cout << "hi\n";

  // Tests
  { using std::cout;
  fProj.Add(geo::Point_t(100, 100, 100), 100, 30);
  cout << fProj.Size() << "\n";

  fProj.ProjectToWires();
  cout << fProj.Size() << " - " << fProj.GetNumInvalidProjections() << "\n";

  for (int i = 0; i < fProj.Size(); i++) {
    for (auto rID : fProj.ActiveROPIDs() ) {
      cout << "rID = " << rID << "\n";
      cout << "ch = " << fProj.GetChs(i)[rID] <<
        " - tick = " << fProj.GetTicks(i)[rID] <<
        " - fd drift = " << fProj.GetFDDrifts(i)[rID] <<
        " - wire distance = " << fProj.GetWireDistances(i)[rID] <<
        " - adc = " << fProj.GetAdc(i) <<
        " - nd drift = " << fProj.GetNDDrift(i) << "\n";
    }
  }

  fProj.Clear();
  cout << fProj.Size() << " - " << fProj.GetNumInvalidProjections() << "\n";
  }
}

void extrapolation::NDToFD::beginJob()
{
  // NOTE Maybe this should all go in the constructor instead for running over multiple files.
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
  // Initialise projector
  // fProj(0.0, 0.0, 0.0, &fGeom, &detProp, true);
  fProj = Projections(0.0, 0.0, 0.0, fGeom, true);

  // Load torchscript models
  try {
    fNetworks[geo::kZ] = torch::jit::load(fNetworkPathZ, torch::kCUDA);
    fNetworks[geo::kZ] = torch::jit::optimize_for_inference(fNetworks[geo::kZ]);
    std::cout << "Z plane module loaded from " << fNetworkPathZ << std::endl;
  }
  catch (const c10::Error& err) {
    std::cerr << "error loading the model\n";
    std::cerr << err.what() << "\n";
  }
  try {
    fNetworks[geo::kU] = torch::jit::load(fNetworkPathU, torch::kCUDA);
    fNetworks[geo::kU] = torch::jit::optimize_for_inference(fNetworks[geo::kU]);
    std::cout << "U plane module loaded from " << fNetworkPathU << std::endl;
  }
  catch (const c10::Error& err) {
    std::cerr << "error loading the model\n";
    std::cerr << err.what() << "\n";
  }
  try {
    fNetworks[geo::kV] = torch::jit::load(fNetworkPathV, torch::kCUDA);
    fNetworks[geo::kV] = torch::jit::optimize_for_inference(fNetworks[geo::kV]);
    std::cout << "V plane module loaded from " << fNetworkPathV << std::endl;
  }
  catch (const c10::Error& err) {
    std::cerr << "error loading the model\n";
    std::cerr << err.what() << "\n";
  }
}

void extrapolation::NDToFD::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::NDToFD)

