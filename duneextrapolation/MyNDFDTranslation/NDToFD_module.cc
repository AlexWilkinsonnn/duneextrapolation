////////////////////////////////////////////////////////////////////////
// Class:       NDToFD
// Plugin Type: producer (Unknown Unknown)
// File:        NDToFD_module.cc
//
// Copied on Tue 21 Jun 22 by Alex Wilkinson.
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
  produces<std::vector<raw::RawDigit>>("NDTranslated");
}

void extrapolation::NDToFD::produce(art::Event& e)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  std::cout << "test\n";

}

void extrapolation::NDToFD::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

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

