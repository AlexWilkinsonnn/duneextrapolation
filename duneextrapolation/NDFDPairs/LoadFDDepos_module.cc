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

#include "highfive/H5DataSet.hpp"
#include "highfive/H5File.hpp"
#include "highfive/H5Group.hpp"
#include "highfive/H5Object.hpp"
#include "highfive/H5DataType.hpp"
#include "highfive/H5DataSpace.hpp"

#include <string>
#include <vector>
#include <map>
#include <set>
#include <math.h>

typedef struct depo {
  int eventID;
  int uniqID;
  double x_end;
  double x_start;
  double y_end;
  double y_start;
  double z_end;
  double z_start;
  double x;
  double y;
  double z;
  double t0_end;
  double t0_start;
  double t0;
  double dx;
  double dEdx;
  double dE;
  int outsideNDLAr;
} depo;

HighFive::CompoundType make_depo() {
  return {
    {"eventID", HighFive::AtomicType<int>{}},
    {"uniqID", HighFive::AtomicType<int>{}},
    {"x_end", HighFive::AtomicType<double>{}},
    {"x_start", HighFive::AtomicType<double>{}},
    {"y_end", HighFive::AtomicType<double>{}},
    {"y_start", HighFive::AtomicType<double>{}},
    {"z_end", HighFive::AtomicType<double>{}},
    {"z_start", HighFive::AtomicType<double>{}},
    {"x", HighFive::AtomicType<double>{}},
    {"y", HighFive::AtomicType<double>{}},
    {"z", HighFive::AtomicType<double>{}},
    {"t0_end", HighFive::AtomicType<double>{}},
    {"t0_start", HighFive::AtomicType<double>{}},
    {"t0", HighFive::AtomicType<double>{}},
    {"dx", HighFive::AtomicType<double>{}},
    {"dEdx", HighFive::AtomicType<double>{}},
    {"dE", HighFive::AtomicType<double>{}},
    {"outsideNDLAr", HighFive::AtomicType<int>{}}
  };
}

typedef struct vertex {
  int eventID;
  double x_vert;
  double y_vert;
  double z_vert;
} depoVtx;

HighFive::CompoundType make_vertex() {
  return {
    {"eventID", HighFive::AtomicType<int>{}},
    {"x_vert", HighFive::AtomicType<double>{}},
    {"y_vert", HighFive::AtomicType<double>{}},
    {"z_vert", HighFive::AtomicType<double>{}},
  };
}

HIGHFIVE_REGISTER_TYPE(depo, make_depo)
HIGHFIVE_REGISTER_TYPE(vertex, make_vertex)

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
  HighFive::File* fFile;
  std::map<int, std::vector<depo>> fDepos;
  std::vector<int> fEventIDs;

  std::string fNDFDH5FileLoc;

  bool fInsideNDOnly;
  double fTargetSEDLen;
};

extrapolation::LoadFDDepos::LoadFDDepos(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fNDFDH5FileLoc (p.get<std::string>("NDFDH5FileLoc")),
    fInsideNDOnly (p.get<bool>("InsideNDOnly")),
    fTargetSEDLen  (p.get<double>("TargetSEDLen"))
{
  produces<std::vector<sim::SimEnergyDeposit>>("LArG4DetectorServicevolTPCActive");
  produces<std::vector<sim::SimEnergyDeposit>>("eventID");
}

void extrapolation::LoadFDDepos::produce(art::Event& e)
{
  if (!fEventIDs.size()) {
    throw cet::exception("LoadFDDepos")
      << "Number of events exceed unique eventIDs in " << fNDFDH5FileLoc
      << " - Line " << __LINE__ << " in file " << __FILE__ << "\n";
  }
  int currEventID = fEventIDs.back();
  fEventIDs.pop_back();
  const std::vector<depo> eventDepos = fDepos[currEventID];

  std::cout << "eventID = " << currEventID << ": Loading " << eventDepos.size() << " depos";
  if (fInsideNDOnly){
    int numDropping = 0;
    for (const depo& dep : eventDepos) {
      if (dep.outsideNDLAr == 1) {
        numDropping++;
      }
    }
    std::cout << " (dropping " << numDropping << ")";
  }
  std::cout << "\n";

  // SED vector with single SED that stores the eventID for matching with ND data later
  // eventID is being stored in the trackID member
  auto evNum = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
  geo::Point_t posStart = geo::Point_t(0, 0, 0);
  geo::Point_t posEnd = geo::Point_t(0, 0, 0);
  sim::SimEnergyDeposit ID = sim::SimEnergyDeposit(
    0, 0, 0, 0, posStart, posEnd, 0, 0, currEventID, 0
  );
  evNum->push_back(ID);

  auto SEDs = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
  for (const depo& dep : eventDepos) {
    if (fInsideNDOnly && dep.outsideNDLAr == 1) {
      continue;
    }
    if (!fTargetSEDLen) { // Load as an SED
      const geo::Point_t posStart(dep.x_start, dep.y_start, dep.z_start);
      const geo::Point_t posEnd(dep.x_end, dep.y_end, dep.z_end);
      sim::SimEnergyDeposit SED = sim::SimEnergyDeposit(
        0, 0, 0, dep.dE, posStart, posEnd, dep.t0_start, dep.t0_end, 1, 1
      );
      SEDs->push_back(SED);
    }
    else { // Break into steps and load each as an SED
      const int nSteps = std::ceil(dep.dx / fTargetSEDLen);
      const double deltaX = (dep.x_end - dep.x_start) / (double)nSteps;
      const double deltaY = (dep.y_end - dep.y_start) / (double)nSteps;
      const double deltaZ = (dep.z_end - dep.z_start) / (double)nSteps;
      const double deltaT0 = (dep.t0_end - dep.t0_start) / (double)nSteps;
      for (int iStep = 0; iStep < nSteps; iStep++) {
        const geo::Point_t posStart(
          dep.x_start + iStep * deltaX,
          dep.y_start + iStep * deltaY,
          dep.z_start + iStep * deltaZ
        );
        const geo::Point_t posEnd(
          dep.x_start + (iStep + 1) * deltaX,
          dep.y_start + (iStep + 1) * deltaY,
          dep.z_start + (iStep + 1) * deltaZ
        );
        // Ideally would account for fluctuations in energy between steps, not possible without
        // running edep-sim with small step size which makes larnd-sim behave weird
        sim::SimEnergyDeposit SED = sim::SimEnergyDeposit(
          0, 0, 0,
          dep.dE / (double)nSteps,
          posStart, posEnd,
          dep.t0_start + iStep * deltaT0, dep.t0_start + (iStep + 1) * deltaT0,
          1, 1
        );
        SEDs->push_back(SED);
      }
    }
  }

  e.put(std::move(SEDs), "LArG4DetectorServicevolTPCActive");
  e.put(std::move(evNum), "eventID");
}

void extrapolation::LoadFDDepos::beginJob()
{
  std::cout << "Readng data from " << fNDFDH5FileLoc << "\n";
  fFile = new HighFive::File(fNDFDH5FileLoc, HighFive::File::ReadOnly);

  // Read in data
  std::vector<vertex> readVtx;
  HighFive::DataSet datasetVtx = fFile->getDataSet("fd_vertices");
  datasetVtx.read(readVtx);
  std::vector<depo> readDepos;
  HighFive::DataSet datasetDepos = fFile->getDataSet("fd_deps");
  datasetDepos.read(readDepos);

  // Create maps eventIDs and data
  for (const depo& dep : readDepos) {
    fDepos[dep.eventID].push_back(dep);
  }

  for (const vertex& vtx : readVtx) {
    fEventIDs.push_back(vtx.eventID);
  }
  if (!fEventIDs.size()) {
    throw cet::exception("LoadFDDepos")
      << "No data to read in " << fNDFDH5FileLoc
      << " - Line " << __LINE__ << " in file " << __FILE__ << "\n";
  }
  std::sort(fEventIDs.begin(), fEventIDs.end(), std::greater<>()); // sort desc so to pop from back
  std::cout << fEventIDs.size() << " unique eventIDs in input\n";
}

void extrapolation::LoadFDDepos::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::LoadFDDepos)
