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
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

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
#include "duneextrapolation/MyNDFDTranslation/Types.h"

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
  std::map<geo::View_t, double>                     fTickShifts;
  std::map<geo::View_t, std::vector<double>>        fInputScaleFactors;
  std::map<geo::View_t, double>                     fOutputScaleFactors;

  Projections fProj;

  // fhicl params
  std::string         fNDPacketsLabel;
  std::string         fNetworkPathZ;
  std::string         fNetworkPathU;
  std::string         fNetworkPathV;
  int                 fPixelMapMode;
  double              fTickShiftZ;
  double              fTickShiftU;
  double              fTickShiftV;
  std::vector<double> fInputScaleFactorsZ;
  std::vector<double> fInputScaleFactorsU;
  std::vector<double> fInputScaleFactorsV;
  double              fOutputScaleFactorZ;
  double              fOutputScaleFactorU;
  double              fOutputScaleFactorV;
  std::vector<int>    fDoubleColZWires;
  bool                fSavePixelMaps;
};

extrapolation::NDToFD::NDToFD(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fNDPacketsLabel     (p.get<std::string>("NDPacketsLabel")),
    fNetworkPathZ       (p.get<std::string>("NetworkPathZ")),
    fNetworkPathU       (p.get<std::string>("NetworkPathU")),
    fNetworkPathV       (p.get<std::string>("NetworkPathV")),
    fPixelMapMode       (p.get<int>("PixelMapMode")),
    fTickShiftZ         (p.get<double>("TickShiftZ")),
    fTickShiftU         (p.get<double>("TickShiftU")),
    fTickShiftV         (p.get<double>("TickShiftV")),
    fInputScaleFactorsZ (p.get<std::vector<double>>("InputScaleFactorsZ")),
    fInputScaleFactorsU (p.get<std::vector<double>>("InputScaleFactorsU")),
    fInputScaleFactorsV (p.get<std::vector<double>>("InputScaleFactorsV")),
    fOutputScaleFactorZ (p.get<double>("OutputScaleFactorZ")),
    fOutputScaleFactorU (p.get<double>("OutputScaleFactorU")),
    fOutputScaleFactorV (p.get<double>("OutputScaleFactorV")),
    fDoubleColZWires    (p.get<std::vector<int>>("DoubleColZWires")),
    fSavePixelMaps      (p.get<bool>("SavePixelMaps"))
{
  consumes<std::vector<sim::SimEnergyDeposit>>(fNDPacketsLabel);

  produces<std::vector<raw::RawDigit>>("NDTranslated");

  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  // Initialise projector
  if (fPixelMapMode == 1) {
    fProj = Projections(fTickShiftZ, fTickShiftU, fTickShiftV, fGeom, true);
  }
  else {
    throw std::runtime_error(std::string("PixelMapMode=") +
      std::to_string(fPixelMapMode) + std::string(" not implemented"));
  }

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

  // Fill maps with fcl params
  fTickShifts[geo::kZ] = fTickShiftZ;
  fTickShifts[geo::kU] = fTickShiftU;
  fTickShifts[geo::kV] = fTickShiftV;
  fInputScaleFactors[geo::kZ] = fInputScaleFactorsZ;
  fInputScaleFactors[geo::kU] = fInputScaleFactorsU;
  fInputScaleFactors[geo::kV] = fInputScaleFactorsV;
  fOutputScaleFactors[geo::kZ] = fOutputScaleFactorZ;
  fOutputScaleFactors[geo::kU] = fOutputScaleFactorU;
  fOutputScaleFactors[geo::kV] = fOutputScaleFactorV;
}

void extrapolation::NDToFD::produce(art::Event& e)
{
  fProj.Clear();

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
  fProj.SetDetProp(&detProp);

  const auto sedsNDPackets = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fNDPacketsLabel);
  auto digs = std::make_unique<std::vector<raw::RawDigit>>();

  for (const sim::SimEnergyDeposit& sed : *sedsNDPackets) {
    fProj.Add(sed);
  }

  fProj.ProjectToWires();

  // Pixel Map mode 1 is the only one at the moment, refactor stuff into functions later
  for (readout::ROPID rID : fProj.ActiveROPIDs()) {
    // Construct ROPID tensor
    torch::Tensor NDTensor = fGeom->View(rID) == geo::kZ ?
      torch::zeros(
        {1, 5, 480, 4492}, torch::dtype(torch::kFloat32).device(torch::kCPU).requires_grad(false)) :
      torch::zeros(
        {1, 6, 800, 4492}, torch::dtype(torch::kFloat32).device(torch::kCPU).requires_grad(false));
    auto NDTensorAccess = NDTensor.accessor<float, 4>();

    std::map<std::pair<float, float>, std::pair<int, std::vector<int>>> pixelTriggers;

    // Fill with data from projetions
    for (ProjectionData proj : fProj.GetProjectionData(rID)) {
      int& ch = proj.localCh;
      int& tick = proj.tick;
      int adc = proj.packet.adc * fInputScaleFactors[fGeom->View(rID)][0];
      // Drift effect goes like sqrt of drift distance. Also doing weighted average by adc.
      double NDDrift = std::sqrt(proj.packet.NDDrift) * fInputScaleFactors[fGeom->View(rID)][1] * adc;
      double FDDrift = std::sqrt(proj.FDDrift) * fInputScaleFactors[fGeom->View(rID)][2] * adc;

      NDTensorAccess[0][0][ch][tick] += adc;
      NDTensorAccess[0][1][ch][tick] += NDDrift;
      NDTensorAccess[0][2][ch][tick] += FDDrift;
      NDTensorAccess[0][3][ch][tick] += fInputScaleFactors[fGeom->View(rID)][3]; // Num stacked

      // Record data on pixel triggers and fill in wire distance channel
      if (fGeom->View(rID) != geo::kZ) {
        float z = std::round((float)proj.packet.pos.Z() * 10.0) / 10.0;
        float y = std::round((float)proj.packet.pos.Y() * 10.0) / 10.0;
        // std::cout << proj.packet.pos.Z() << " - " << z << "\n " << proj.packet.pos.Y() << " - " << y << "\n";

        if (!pixelTriggers.count(std::make_pair(z, y))) {
          pixelTriggers[std::make_pair(z, y)] = std::make_pair(ch, std::vector<int> { tick });
        }
        else {
          pixelTriggers[std::make_pair(z, y)].second.push_back(tick);
        }

        double wireDistance = proj.wireDistance * fInputScaleFactors[fGeom->View(rID)][5] * adc;

        NDTensorAccess[0][5][ch][tick] += wireDistance;
      }
      // Fill double pixel column channels
      else {
        for (int doubleColCh : fDoubleColZWires) {
          NDTensor.index_put_({0, 4, doubleColCh, torch::indexing::Slice()},
            fInputScaleFactors[fGeom->View(rID)][4]);
        }
      }
    }

    // Last step in the adc weighted average
    for (ProjectionData proj : fProj.GetProjectionData(rID)) {
      if (NDTensorAccess[0][0][proj.localCh][proj.tick] != 0) {
        NDTensorAccess[0][1][proj.localCh][proj.tick] /=
          NDTensorAccess[0][0][proj.localCh][proj.tick];
        NDTensorAccess[0][2][proj.localCh][proj.tick] /=
          NDTensorAccess[0][0][proj.localCh][proj.tick];

        if (fGeom->View(rID) != geo::kZ) {
          NDTensorAccess[0][5][proj.localCh][proj.tick] /=
            NDTensorAccess[0][0][proj.localCh][proj.tick];
        }
      }
    }

    // Use pixel trigger data to identify first triggers and fill in tensor with first trigger data.
    if (fGeom->View(rID) != geo::kZ) {
      for (std::pair<std::pair<int, int>, std::pair<int, std::vector<int>>> pixelTrigger : pixelTriggers ) {
        std::sort(pixelTrigger.second.second.begin(), pixelTrigger.second.second.end()); // Sort in place

        // Identify first tiggers in self-triggering cycle
        std::vector<int> firstTriggerTicks;
        for (int i = 0; i < (int)pixelTrigger.second.second.size(); i++) {
          if (i == 0 || pixelTrigger.second.second[i] - firstTriggerTicks[i - 1] > 15) {
            firstTriggerTicks.push_back(pixelTrigger.second.second[i]);
          }
        }

        // Fill tensor with first triggers.
        for (int tick : firstTriggerTicks) {
          NDTensorAccess[0][4][pixelTrigger.second.first][tick] +=
            fInputScaleFactors[fGeom->View(rID)][4];
        }
      }
    }

    // Do inference
    torch::NoGradGuard no_grad_guard;

    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(NDTensor.to(torch::kCUDA));

    torch::Tensor FDTranslatedTensor = fNetworks[fGeom->View(rID)].forward(inputs).toTensor().detach().to(torch::kCPU);

    FDTranslatedTensor = FDTranslatedTensor[0][0] / fOutputScaleFactors[fGeom->View(rID)];
    FDTranslatedTensor = FDTranslatedTensor.to(torch::kShort);

    // Write out to raw digits
    auto FDTranslatedTensorAccess = FDTranslatedTensor.accessor<short, 2>();
    for (int i = 0; i < (int)fGeom->Nchannels(rID); i++) {
      raw::RawDigit::ADCvector_t adcVec(4492);
      for (int j = 0; j < 4492; j++) {
        adcVec[j] = FDTranslatedTensorAccess[i][j];
      }

      raw::RawDigit dig((raw::ChannelID_t)i + fGeom->FirstChannelInROP(rID), adcVec.size(), adcVec);
      dig.SetPedestal(0);
      digs->push_back(dig);
    }
  }

  // Fill remaining channels with zero vectors
  std::vector<readout::ROPID> rIDsWritten = fProj.ActiveROPIDs();
  for (readout::ROPID rID : fGeom->IterateROPIDs()) {
    if (std::find(rIDsWritten.begin(), rIDsWritten.end(), rID) == rIDsWritten.end()) {
      raw::ChannelID_t firstCh = fGeom->FirstChannelInROP(rID);
      for (raw::ChannelID_t ch = firstCh; ch < firstCh + fGeom->Nchannels(rID); ch++) {
        raw::RawDigit::ADCvector_t adcVec(4492, 0);

        raw::RawDigit dig(ch, adcVec.size(), adcVec);
        dig.SetPedestal(0);
        digs->push_back(dig);
      }
    }
  }

  e.put(std::move(digs), "NDTranslated");
}

void extrapolation::NDToFD::beginJob()
{
}

void extrapolation::NDToFD::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::NDToFD)

