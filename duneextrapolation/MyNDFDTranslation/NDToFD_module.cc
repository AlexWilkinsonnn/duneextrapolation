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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// #include <memory>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include <torch/torch.h>
#include <torch/script.h>

// All ROOT stuff needs to go after torch stuff because of the ClassDef macro ROOT defines.
#include "TH2D.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

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
  void fillNDTensor(torch::Tensor& NDTensor, readout::ROPID rID);

  void addNDHits(art::Event& e);

  void writeTensorToTH2(torch::Tensor tensor, std::string name, const double scaleFactor = 1.0);
  void writeDigitsToTH2(art::Event& e, std::vector<readout::ROPID> rIDs);

  void applySignalMask(torch::Tensor NDAdcTensor, torch::Tensor& FDTensor, const int maxChShift,
    const int maxTickShift);

  const geo::GeometryCore* fGeom;

  std::map<geo::View_t, torch::jit::script::Module> fNetworks;
  std::map<geo::View_t, double>                     fTickShifts;
  std::map<geo::View_t, std::vector<double>>        fInputScaleFactors;
  std::map<geo::View_t, double>                     fOutputScaleFactors;
  std::map<geo::View_t, std::map<ChannelType, int>> fExtraChannels;

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
  bool                fSaveTrueDigits;
  std::string         fTrueDigitsLabel;
  bool                fMakeNDHits;
  double              fNDHitsScaleFactor;
  bool                fApplySignalMask;
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
    fSavePixelMaps      (p.get<bool>("SavePixelMaps")),
    fSaveTrueDigits     (p.get<bool>("SaveTrueDigits")),
    fTrueDigitsLabel    (p.get<std::string>("TrueDigitsLabel")),
    fMakeNDHits         (p.get<bool>("MakeNDHits")),
    fNDHitsScaleFactor  (p.get<double>("NDHitsScaleFactor")),
    fApplySignalMask    (p.get<bool>("ApplySignalMask"))
{
  consumes<std::vector<sim::SimEnergyDeposit>>(fNDPacketsLabel);
  if (fSaveTrueDigits) {
    consumes<std::vector<raw::RawDigit>>(fTrueDigitsLabel);
  }

  produces<std::vector<raw::RawDigit>>("NDTranslated");
  if (fMakeNDHits) {
    produces<std::vector<recob::Hit>>("NDPackets");
  }

  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  // Initialise projector
  if (fPixelMapMode == 1) {
    fProj = Projections(fTickShiftZ, fTickShiftU, fTickShiftV, fGeom, true);
    fExtraChannels[geo::kZ] = { { kDoubleColZWires, 4 } };
    fExtraChannels[geo::kU] = { { kFirstTriggers, 4 }, { kWireDistance, 5 } };
    fExtraChannels[geo::kV] = { { kFirstTriggers, 4 }, { kWireDistance, 5 } };
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
  auto digs = std::make_unique<std::vector<raw::RawDigit>>();

  fProj.Clear();

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
  fProj.SetDetProp(&detProp);

  const auto sedsNDPackets = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fNDPacketsLabel);
  for (const sim::SimEnergyDeposit& sed : *sedsNDPackets) {
    fProj.Add(sed);
  }

  fProj.ProjectToWires();

  for (readout::ROPID rID : fProj.ActiveROPIDs()) {
    const geo::View_t view = fGeom->View(rID);
    const int nChs = (int)fGeom->Nchannels(rID);

    torch::Tensor NDInput = torch::zeros({1, 4 + (int)fExtraChannels[view].size(), nChs, 4492},
      torch::dtype(torch::kFloat32).device(torch::kCPU).requires_grad(false));

    // Fill ND tensor with projection information
    fillNDTensor(NDInput, rID);

    // Do inference
    torch::NoGradGuard no_grad_guard;

    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(NDInput.to(torch::kCUDA));

    torch::Tensor FDOutput = fNetworks[view].forward(inputs).toTensor().detach().to(torch::kCPU);

    FDOutput = FDOutput[0][0] / fOutputScaleFactors[view];

    if (fApplySignalMask) {
      applySignalMask(NDInput.to(torch::kCPU)[0][0] / fInputScaleFactors[view][0], FDOutput, 2, 25);
    }

    FDOutput = FDOutput.to(torch::kShort);

    // Write out to raw digits
    auto FDOutputAccess = FDOutput.accessor<short, 2>();
    for (int i = 0; i < nChs; i++) {
      raw::RawDigit::ADCvector_t adcVec(4492);
      for (int j = 0; j < 4492; j++) {
        adcVec[j] = FDOutputAccess[i][j];
      }

      raw::RawDigit dig((raw::ChannelID_t)i + fGeom->FirstChannelInROP(rID), adcVec.size(), adcVec);
      dig.SetPedestal(0);
      digs->push_back(dig);
    }

    if (fSavePixelMaps) {
      std::string FDName = "pm_fd_rop" + std::to_string(rID.ROP) + "_tpcset" +
        std::to_string(rID.TPCset) + "_view" + std::to_string(fGeom->View(rID)) + "_e" +
        std::to_string(e.event());
      writeTensorToTH2(FDOutput.to(torch::kFloat), FDName);

      for (int iDepth = 0; iDepth < NDInput.sizes()[1]; iDepth++) {
        std::string NDName = "pm_nd_rop" + std::to_string(rID.ROP) + "_tpcset" +
          std::to_string(rID.TPCset) + "_view" + std::to_string(fGeom->View(rID)) + "_e" +
          std::to_string(e.event()) + "_ch" + std::to_string(iDepth);
        writeTensorToTH2(NDInput[0][iDepth], NDName, fInputScaleFactors[view][iDepth]);
      }
    }
  } // for (readout::ROPID rID : fProj.ActiveROPIDs())

  if (fMakeNDHits) {
    addNDHits(e);
  }

  if (fSaveTrueDigits) {
    writeDigitsToTH2(e, fProj.ActiveROPIDs());
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

// Fill tensor with the default adc, drifts, stacked then request extra channels using a map of
// channel number to type of channel.
void extrapolation::NDToFD::fillNDTensor(torch::Tensor& NDTensor, readout::ROPID rID)
{
  std::map<std::pair<float, float>, std::pair<int, std::vector<int>>> pixelTriggers;
  const geo::View_t view = fGeom->View(rID);
  std::map<ChannelType, int> extraChannels = fExtraChannels[view];

  auto NDTensorAccess = NDTensor.accessor<float, 4>();

  for (ProjectionData proj : fProj.GetProjectionData(rID)) {
    const int ch = proj.localCh, tick = proj.tick;
    double adc = (double)proj.packet.adc * fInputScaleFactors[view][0];
    // Drift effect goes like sqrt of drift distance. Also doing weighted average by adc.
    double NDDrift = std::sqrt(proj.packet.NDDrift) * fInputScaleFactors[view][1] * adc;
    double FDDrift = std::sqrt(proj.FDDrift) * fInputScaleFactors[view][2] * adc;

    // Channels common to all ND tensors
    NDTensorAccess[0][0][ch][tick] += adc;
    NDTensorAccess[0][1][ch][tick] += NDDrift;
    NDTensorAccess[0][2][ch][tick] += FDDrift;
    NDTensorAccess[0][3][ch][tick] += fInputScaleFactors[view][3]; // Num stacked

    if (extraChannels.count(kFirstTriggers)) {
      float z = std::round((float)proj.packet.pos.Z() * 10.0) / 10.0;
      float y = std::round((float)proj.packet.pos.Y() * 10.0) / 10.0;

      if (!pixelTriggers.count(std::make_pair(z, y))) {
        pixelTriggers[std::make_pair(z, y)] = std::make_pair(ch, std::vector<int> { tick });
      }
      else {
        pixelTriggers[std::make_pair(z, y)].second.push_back(tick);
      }
    }

    if (extraChannels.count(kWireDistance)) {
      double wireDistance = proj.wireDistance *
        fInputScaleFactors[view][extraChannels[kWireDistance]] * adc;

      NDTensorAccess[0][extraChannels[kWireDistance]][ch][tick] += wireDistance;
    }

    if (extraChannels.count(kDoubleColZWires)) {
      for (int doubleColCh : fDoubleColZWires) {
        NDTensor.index_put_({0, 4, doubleColCh, torch::indexing::Slice()},
          fInputScaleFactors[view][extraChannels[kDoubleColZWires]]);
      }
    }
  }

  // Last step in the adc weighted average
  for (std::pair<int, int> chTick : fProj.ActiveChTicks(rID)) {
    const int ch = chTick.first, tick = chTick.second;

    if (NDTensorAccess[0][0][ch][tick] != 0) {
      NDTensorAccess[0][1][ch][tick] /=
        NDTensorAccess[0][0][ch][tick];
      NDTensorAccess[0][2][ch][tick] /=
        NDTensorAccess[0][0][ch][tick];

      if (extraChannels.count(kWireDistance)) {
        NDTensorAccess[0][extraChannels[kWireDistance]][ch][tick] /=
          NDTensorAccess[0][0][ch][tick];
      }
    }
  }

  if (extraChannels.count(kFirstTriggers)) {
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
        NDTensorAccess[0][extraChannels[kFirstTriggers]][pixelTrigger.second.first][tick] +=
          fInputScaleFactors[view][extraChannels[kFirstTriggers]];
      }
    }
  }
}

// Make hits directly out of ND packets and put them on event
void extrapolation::NDToFD::addNDHits(art::Event& e)
{
  auto hits = std::make_unique<std::vector<recob::Hit>>();

  for (readout::ROPID rID : fProj.ActiveROPIDs()) {
    const geo::View_t view = fGeom->View(rID);

    for (ProjectionData proj : fProj.GetProjectionData(rID)) {
      raw::ChannelID_t ch = (raw::ChannelID_t)proj.localCh + fGeom->FirstChannelInROP(rID);

      float peakTime = (float)proj.tick;
      raw::TDCtick_t startTick = (raw::TDCtick_t)(peakTime);
      raw::TDCtick_t endTick = (raw::TDCtick_t)(peakTime + 5.0);

      // empirical scaling to be ~ FD scale
      float integral = (float)((double)(proj.packet.adc * 16) * fNDHitsScaleFactor);
      float peakAmplitude = integral / 5.0;
      float summedADC = (float)proj.packet.adc * 16;

      float rms = 3.0; // fairly arbitrary choice

      // Gauss hit finder just takes the first wire, I guess the disaambiguation modules will get
      // correct wireID later
      // (https://internal.dunescience.org/doxygen/GausHitFinder__module_8cc_source.html)
      geo::WireID wID = fGeom->ChannelToWire(ch)[0];

      recob::Hit hit(ch, startTick, endTick, peakTime, float(-1.0), rms, peakAmplitude, float(-1.0),
        summedADC, integral, float(-1.0), short(0), short(-1), float(0.0), int(-1), view,
        geo::SigType_t(geo::kMysteryType), wID);

      hits->push_back(hit);
    }
  }

  e.put(std::move(hits), "NDPackets");
}


void extrapolation::NDToFD::writeTensorToTH2(torch::Tensor tensor, std::string name,
  const double scaleFactor /* = 1.0 */)
{
  art::ServiceHandle<art::TFileService> tfs;
  auto tensorAccess = tensor.accessor<float, 2>();

  int nChs = tensor.sizes()[0];

  TH2D* h = tfs->make<TH2D>(name.c_str(), name.c_str(), nChs, 0, nChs, 4492, 0 , 4492);

  for(int iCh = 0; iCh < nChs; iCh++) {
    for (int iTick = 0; iTick < 4492; iTick++) {
      h->SetBinContent(iCh + 1, iTick + 1, tensorAccess[iCh][iTick] / scaleFactor);
    }
  }

  h->Write();

  delete h;
}

void extrapolation::NDToFD::writeDigitsToTH2(art::Event& e, std::vector<readout::ROPID> rIDs)
{
  art::ServiceHandle<art::TFileService> tfs;
  const auto digs = e.getValidHandle<std::vector<raw::RawDigit>>(fTrueDigitsLabel);

  std::map<readout::ROPID, TH2D*> hs;
  for (readout::ROPID rID : rIDs) {
    std::string name = "pm_fdtrue_rop" + std::to_string(rID.ROP) + "_tpcset" +
      std::to_string(rID.TPCset) + "_view" + std::to_string(fGeom->View(rID)) + "_e" +
      std::to_string(e.event());
    hs[rID] = tfs->make<TH2D>(name.c_str(), name.c_str(), fGeom->Nchannels(rID), 0,
      fGeom->Nchannels(rID), 4492, 0, 4492);
  }

  for (raw::RawDigit dig : *digs) {
    readout::ROPID rID = fGeom->ChannelToROP(dig.Channel());

    if (hs.count(rID)) {
      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

      // WC might still have a 6000 tick readout window but just ignoring everything after
      // 4492
      for (int iTick = 0; iTick < 4492; iTick++) {
        const short adc = adcs[iTick] ? short(adcs[iTick]) - dig.GetPedestal() : 0;

        hs[rID]->SetBinContent(dig.Channel() - fGeom->FirstChannelInROP(rID) + 1, iTick + 1, adc);
      }
    }
  }

  for (auto pair : hs) {
    pair.second->Write();
    delete pair.second;
  }
}

void extrapolation::NDToFD::applySignalMask(torch::Tensor NDAdcTensor, torch::Tensor& FDTensor,
  const int maxChShift, const int maxTickShift)
{
  // clone to prevent overlapping memory erros in later operations
  torch::Tensor mask = NDAdcTensor.clone();

  // Make signal mask by smearing ND tensor
  { using namespace torch::indexing;
  for (int tickShift = 1; tickShift <= maxTickShift; tickShift++) {
    mask.index({Slice(), Slice(tickShift, None)}) +=
      NDAdcTensor.index({Slice(), Slice(None, -tickShift)});
    mask.index({Slice(), Slice(None, -tickShift)}) +=
      NDAdcTensor.index({Slice(), Slice(tickShift, None)});
  }

  for (int chShift = 1; chShift <= maxChShift; chShift++) {
    mask.index({Slice(chShift, None), Slice()}) +=
      NDAdcTensor.index({Slice(None, -chShift), Slice()});
    mask.index({Slice(None, -chShift), Slice()}) +=
      NDAdcTensor.index({Slice(chShift, None), Slice()});
  }
  } // using namespace torch::indexing;

  mask = mask.to(torch::kBool).to(torch::kFloat32);

  FDTensor.mul_(mask);
}

DEFINE_ART_MODULE(extrapolation::NDToFD)

