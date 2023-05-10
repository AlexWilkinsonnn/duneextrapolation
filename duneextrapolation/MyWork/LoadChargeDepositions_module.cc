////////////////////////////////////////////////////////////////////////
// Class:       LoadChargeDepositions
// Plugin Type: producer (Unknown Unknown)
// File:        LoadChargeDepositions.cc
//
// Generated at May 10 2023 by Alexander Wilkinson.
//
// Read charge deposition information from another source into a
// SimEnergyDeposit product. Expects the fhicl to point to a root with
// a TTree called nd_depos that has the format for each event:
// nd_depos : vector<double array[12]>
//   array[0] - track id
//   array[1] - pdg code
//   array[2] - start X (cm)
//   array[3] - end X (cm)
//   array[4] - start Y (cm)
//   array[5] - end Y (cm)
//   array[6] - start Z (cm)
//   array[7] - end Z (cm)
//   array[8] - start T (us)
//   array[9] - end T (us)
//   array[10] - num electrons
//   array[11] - dE
// vertex : double array[4]
//   array[0] - X (cm)
//   array[0] - Y (cm)
//   array[0] - Z (cm)
//   array[0] - T (us)
// eventID : int
// IMPORTANT Above coordinates are in the FD convention where X is the drift
// direction and Z is beam direction. ND convention is to swap these.
//
// NOTE Currently expect the ionistaion calculation to be done already
// to get number of electrons. Might be better to have IonAndScint do
// this instead and just load in dE
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
  class LoadChargeDepositions;
}

class extrapolation::LoadChargeDepositions : public art::EDProducer {
public:
  explicit LoadChargeDepositions(fhicl::ParameterSet const& p);

  LoadChargeDepositions(LoadChargeDepositions const&) = delete;
  LoadChargeDepositions(LoadChargeDepositions&&) = delete;
  LoadChargeDepositions& operator=(LoadChargeDepositions const&) = delete;
  LoadChargeDepositions& operator=(LoadChargeDepositions&&) = delete;

  void produce(art::Event& e) override;

  void beginJob() override;
  void endJob() override;

  void reset();

private:
  const geo::GeometryCore* fGeom;

  // Reading input ND tree
  int                               fNEntries;
  int                               fEntry;
  TTree*                            fTreeDepos;
  std::vector<std::vector<double>>* fDepos;
  std::vector<double>*              fVertex;

  // fhicl params
  std::string fDepoDataLoc;
  double      fXShift;
  double      fYShift;
  double      fZShift;
  double      fZCutLow;
  double      fZCutHigh;

   // Other members
  int            fEventNumber;
};

extrapolation::LoadChargeDepositions::LoadChargeDepositions(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fDepoDataLoc (p.get<std::string>("DepoDataLoc")),
    fXShift      (p.get<double>("XShift")),
    fYShift      (p.get<double>("YShift")),
    fZShift      (p.get<double>("ZShift")),
    fZCutLow     (p.get<double>("ZCutLow")),
    fZCutHigh    (p.get<double>("ZCutHigh"))
{
  produces<std::vector<sim::SimEnergyDeposit>>();
  produces<std::vector<sim::SimEnergyDeposit>>("EventNumber");
}

void extrapolation::LoadChargeDepositions::produce(art::Event& e)
{
  if (fEntry >= fNEntries) {
    std::cout << "Gone beyond number of entries in tree (" << fNEntries << ")\n";
    return;
  }

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  fTreeDepos->GetEntry(fEntry);

  // Make and add the ND depos
  auto SEDs = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
  auto evNum = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

  // For moving events from ND LAr, using YShift=200 the rest 0
  for (const std::vector<double>& depo : *fDepos) {
    int trackID = (int)depo[0];
    int pdg = (int)depo[1];
    double xMin = depo[2] + fXShift;
    double xMax = depo[3] + fXShift;
    double yMin = depo[4] + fYShift;
    double yMax = depo[5] + fYShift;
    double zMin = depo[6] + fZShift;
    double zMax = depo[7] + fZShift;
    double tMin = depo[8];
    double tMax = depo[9];
    int electrons = (int)depo[10];
    double dE = depo[11];

    if ((zMin < fZCutLow || zMax < fZCutLow) || (zMin > fZCutHigh || zMax > fZCutHigh)) {
      continue;
    }

    geo::Point_t posStart = geo::Point_t(xMin, yMin, zMin);
    geo::Point_t posEnd = geo::Point_t(xMax, yMax, zMax);
    
    // if (
    //   (xMin < -363.376 || xMin > 363.376) ||
    //   (yMin < 0 || yMin > 607.829) ||
    //   (zMin < -0.87625 || zMin > 1393.46)
    // ) {
    //   std::cout << "!!!!!!!!!!" << xMin << "," << yMin << "," << zMin << "!!!!!!!!!!\n";
    // }

    sim::SimEnergyDeposit SED = sim::SimEnergyDeposit(
      0, electrons, 0, dE, posStart, posEnd, tMin, tMax, trackID, pdg
    );
    SEDs->push_back(SED);
  }

  // SED that stores an ID for this event
  geo::Point_t posStart = geo::Point_t(0,0,0);
  geo::Point_t posEnd = geo::Point_t(0,0,0);
  sim::SimEnergyDeposit ID = sim::SimEnergyDeposit(
    0, 0, 0, 0, posStart, posEnd, 0, 0, fEventNumber
  );
  evNum->push_back(ID);

  e.put(std::move(SEDs));
  e.put(std::move(evNum), "EventNumber");

  fEventNumber++;
  fEntry++;
}

void extrapolation::LoadChargeDepositions::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  std::cout << "Reading file from " << fDepoDataLoc << "\n";
  fEntry = 0; fDepos = nullptr; fVertex = nullptr;
  TFile* fileDepos = new TFile(fDepoDataLoc.c_str());
  fTreeDepos = (TTree*)fileDepos->Get("nd_depos");
  fTreeDepos->SetBranchAddress("nd_depos", &fDepos);
  fTreeDepos->SetBranchAddress("vertex", &fVertex);

  fEventNumber = 0;
  fNEntries = fTreeDepos->GetEntries();
  std::cout << "File has " << fNEntries << " entries\n";
}

void extrapolation::LoadChargeDepositions::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::LoadChargeDepositions)

