////////////////////////////////////////////////////////////////////////
// Class:       ExportCVNMaps
// Plugin Type: analyzer (Unknown Unknown)
// File:        ExportCVNMaps_module.cc
//
// Wed Apr 13 2022 Alex Wilkinson
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

#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

#include "dunereco/CVN/func/PixelMap.h"
#include "dunereco/CVN/func/Result.h"

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "TH2F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace extrapolation {
  class ExportCVNMaps;
}


class extrapolation::ExportCVNMaps : public art::EDAnalyzer {
public:
  explicit ExportCVNMaps(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ExportCVNMaps(ExportCVNMaps const&) = delete;
  ExportCVNMaps(ExportCVNMaps&&) = delete;
  ExportCVNMaps& operator=(ExportCVNMaps const&) = delete;
  ExportCVNMaps& operator=(ExportCVNMaps&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My functions.
  void reset();

private:
  const geo::GeometryCore* fGeom;

  TTree* fTreePixelMaps;
  int    fEventNum;
  TH2F   fhTruePixelMap;
  TH2F   fhNDPixelMap;
  TH2F   fhNetworkPixelMap;
  // Only interested in flavour scores for now
  float  fTrueNCScore;
  float  fNDNCScore;
  float  fNetworkNCScore;
  float  fTrueNumuScore;
  float  fNDNumuScore;
  float  fNetworkNumuScore;
  float  fTrueNueScore;
  float  fNDNueScore;
  float  fNetworkNueScore;

  std::vector<int> fEventNums;

  std::string fTruePixelMapLabel;
  std::string fNDPixelMapLabel;
  std::string fNetworkPixelMapLabel;
  std::string fTrueCVNResultsLabel;
  std::string fNDCVNResultsLabel;
  std::string fNetworkCVNResultsLabel;
};


extrapolation::ExportCVNMaps::ExportCVNMaps(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fEventNums              (p.get<std::vector<int>> ("EventNumbers")),
    fTruePixelMapLabel      (p.get<std::string>      ("TruePixelMapLabel")),
    fNDPixelMapLabel        (p.get<std::string>      ("NDPixelMapLabel")),
    fNetworkPixelMapLabel   (p.get<std::string>      ("NetworkPixelMapLabel")),
    fTrueCVNResultsLabel    (p.get<std::string>      ("TrueCVNResultsLabel")),
    fNDCVNResultsLabel      (p.get<std::string>      ("NDCVNResultsLabel")),
    fNetworkCVNResultsLabel (p.get<std::string>      ("NetworkCVNResultsLabel"))
{
  consumes<std::vector<cvn::PixelMap>>(fTruePixelMapLabel);
  consumes<std::vector<cvn::PixelMap>>(fNDPixelMapLabel);
  consumes<std::vector<cvn::PixelMap>>(fNetworkPixelMapLabel);

  consumes<std::vector<cvn::Result>>(fTrueCVNResultsLabel);
  consumes<std::vector<cvn::Result>>(fNetworkCVNResultsLabel);
  consumes<std::vector<cvn::Result>>(fNDCVNResultsLabel);

  art::ServiceHandle<art::TFileService> tfs;

  fTreePixelMaps = tfs->make<TTree>("pixelmaps", "Pixel Maps");
  fTreePixelMaps->Branch("EventNum", &fEventNum, "eventnum/I");
  fTreePixelMaps->Branch("TruePixelMap", &fhTruePixelMap);
  fTreePixelMaps->Branch("NDPixelMap", &fhNDPixelMap);
  fTreePixelMaps->Branch("NetworkPixelMap", &fhNetworkPixelMap);
  fTreePixelMaps->Branch("TrueNumuScore", &fTrueNumuScore, "truenumuscore/F");
  fTreePixelMaps->Branch("NetworkNumuScore", &fNetworkNumuScore, "networknumuscore/F");
  fTreePixelMaps->Branch("NDNumuScore", &fNDNumuScore, "ndnumuscore/F");
  fTreePixelMaps->Branch("TrueNueScore", &fTrueNueScore, "truenuescore/F");
  fTreePixelMaps->Branch("NetworkNueScore", &fNetworkNueScore, "networknuescore/F");
  fTreePixelMaps->Branch("NDNueScore", &fNDNueScore, "ndnuescore/F");  
  fTreePixelMaps->Branch("TrueNCScore", &fTrueNCScore, "truencscore/F");
  fTreePixelMaps->Branch("NetworkNCScore", &fNetworkNCScore, "networkncscore/F");
  fTreePixelMaps->Branch("NDNCScore", &fNDNCScore, "ndncscore/F");
}

void extrapolation::ExportCVNMaps::analyze(art::Event const& e)
{
  const auto truePixelMap = e.getValidHandle<std::vector<cvn::PixelMap>>(fTruePixelMapLabel);
  const auto NDPixelMap = e.getValidHandle<std::vector<cvn::PixelMap>>(fNDPixelMapLabel);
  const auto networkPixelMap = e.getValidHandle<std::vector<cvn::PixelMap>>(fNetworkPixelMapLabel);

  const auto trueCVNResults = e.getValidHandle<std::vector<cvn::Result>>(fTrueCVNResultsLabel);
  const auto networkCVNResults = e.getValidHandle<std::vector<cvn::Result>>(fNetworkCVNResultsLabel);
  const auto NDCVNResults = e.getValidHandle<std::vector<cvn::Result>>(fNDCVNResultsLabel);

  if(std::find(fEventNums.begin(), fEventNums.end(), e.id().event()) != fEventNums.end()) {
    this->reset();

    fEventNum = e.id().event();

    fTrueNumuScore = trueCVNResults->at(0).GetNumuProbability();
    fNetworkNumuScore = networkCVNResults->at(0).GetNumuProbability();
    fNDNumuScore = NDCVNResults->at(0).GetNumuProbability();
    fTrueNueScore = trueCVNResults->at(0).GetNueProbability();
    fNetworkNueScore = networkCVNResults->at(0).GetNueProbability();
    fNDNueScore = NDCVNResults->at(0).GetNueProbability();
    fTrueNCScore = trueCVNResults->at(0).GetNCProbability();
    fNetworkNCScore = networkCVNResults->at(0).GetNCProbability();
    fNDNCScore = NDCVNResults->at(0).GetNCProbability();
    
    fhTruePixelMap = *((TH2F*)truePixelMap->at(0).ToTH2()->Clone());
    fhNDPixelMap = *((TH2F*)NDPixelMap->at(0).ToTH2()->Clone());
    fhNetworkPixelMap = *((TH2F*)networkPixelMap->at(0).ToTH2()->Clone());

    fTreePixelMaps->Fill();
  }
}

void extrapolation::ExportCVNMaps::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
}

void extrapolation::ExportCVNMaps::endJob()
{
}

void extrapolation::ExportCVNMaps::reset()
{
  fEventNum = -999;

  fhTruePixelMap = TH2F();
  fhNDPixelMap = TH2F();
  fhNetworkPixelMap = TH2F();

  fTrueNumuScore = -999.0;
  fNetworkNumuScore = -999.0;
  fNDNumuScore = -999.0;
  fTrueNueScore = -999.0;
  fNetworkNueScore = -999.0;
  fNDNueScore = -999.0;
  fTrueNCScore = -999.0;
  fNetworkNCScore = -999.0;
  fNDNCScore = -999.0;
}

DEFINE_ART_MODULE(extrapolation::ExportCVNMaps)
