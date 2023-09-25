////////////////////////////////////////////////////////////////////////
// Class:       LoadFDDepos
// Plugin Type: producer (Unknown Unknown)
// File:        LoadFDDepos_module.cc
//
// Crated on 25 Sep 23 Alex Wilkinson
// Load FD depos from hdf5 file
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

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "hep_hpc/hdf5/File.hpp"
#include "hep_hpc/hdf5/make_ntuple.hpp"

#include <string>
#include <vector>
#include <map>

namespace extrapolation {
  class LoadFDDepos;
}


class extrapolation::LoadFDDepos : public art::EDProducer {
public:
  explicit LoadFDDepos(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  // The line "hep_hpc::hdf5::File fFile;" somehow causes the destructor to have noexcept(false)
  // then there is a build error "looser throw specifier ..." coming from the destructor being
  // noexcept in ED{Analyzer/Producer}. Defining the destructor as below is the only way I found
  // to remove  build error. duhduhduhhhhh code????
  virtual ~LoadFDDepos() noexcept { };

  // Plugins should not be copied or assigned.
  LoadFDDepos(LoadFDDepos const&) = delete;
  LoadFDDepos(LoadFDDepos&&) = delete;
  LoadFDDepos& operator=(LoadFDDepos const&) = delete;
  LoadFDDepos& operator=(LoadFDDepos&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My function.
  void reset();

private:
  const geo::GeometryCore* fGeom;

  hep_hpc::hdf5::File fFile;

  std::string fNDFDH5FileLoc;
};

extrapolation::LoadFDDepos::LoadFDDepos(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fNDFDH5FileLoc (p.get<std::string>("NDFDH5FileLoc"))
{
  produces<std::vector<sim::SimEnergyDeposit>>("LArG4DetectorServicevolTPCActive");
}

void extrapolation::LoadFDDepos::produce(art::Event& e)
{
  std::cout << fFile["fd_deps"]["eventID"][0] << "\n";

  auto SEDs = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
  e.put(std::move(SEDs), "LArG4DetectorServicevolTPCActive");
}

void extrapolation::LoadFDDepos::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  fFile = hep_hpc::hdf5::File(fNDFDH5FileLoc, H5F_ACC_RDONLY);
}

void extrapolation::LoadFDDepos::endJob()
{
  fFile.close();
}

DEFINE_ART_MODULE(extrapolation::LoadFDDepos)
