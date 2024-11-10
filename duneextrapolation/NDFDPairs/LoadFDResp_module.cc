////////////////////////////////////////////////////////////////////////
// Class:       LoadFDResp
// Plugin Type: producer (Unknown Unknown)
// File:        LoadFDResp_module.cc
//
// Crated on 10 Nov 24 Alex Wilkinson
// Load FD wire plane response from hdf5 file
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
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/raw.h"

#include "highfive/H5DataSet.hpp"
#include "highfive/H5File.hpp"
#include "highfive/H5Group.hpp"
#include "highfive/H5Object.hpp"
#include "highfive/H5DataType.hpp"
#include "highfive/H5DataSpace.hpp"

#include <string>
#include <vector>

namespace extrapolation {
  class LoadFDResp;
}

class extrapolation::LoadFDResp : public art::EDProducer {
public:
  explicit LoadFDResp(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  LoadFDResp(LoadFDResp const&) = delete;
  LoadFDResp(LoadFDResp&&) = delete;
  LoadFDResp& operator=(LoadFDResp const&) = delete;
  LoadFDResp& operator=(LoadFDResp&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My function.
  void reset();

private:
  // Members
  const geo::GeometryCore* fGeom;

  HighFive::File* fFile;
  std::vector<int> fEventIDs;

  // Set from fcl
  std::string fNDFDH5FileLoc;
  unsigned int fDigTickWindow;
};

extrapolation::LoadFDResp::LoadFDResp(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fNDFDH5FileLoc (p.get<std::string>("NDFDH5FileLoc")),
    fDigTickWindow (p.get<unsigned int>("DigTickWindow"))
{
  produces<std::vector<sim::SimEnergyDeposit>>("eventID");
  produces<std::vector<raw::RawDigit>>("NDTranslated");
}

void extrapolation::LoadFDResp::produce(art::Event& e)
{
  if (!fEventIDs.size()) {
    throw cet::exception("LoadFDResp")
      << "Number of events exceed unique eventIDs in " << fNDFDH5FileLoc
      << " - Line " << __LINE__ << " in file " << __FILE__ << "\n";
  }
  const int currEventID = fEventIDs.back();
  fEventIDs.pop_back();

  // SED vector with single SED that stores the eventID for matching with ND data later
  // eventID is being stored in the trackID member
  auto evNum = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
  geo::Point_t posStart = geo::Point_t(0, 0, 0);
  geo::Point_t posEnd = geo::Point_t(0, 0, 0);
  sim::SimEnergyDeposit ID = sim::SimEnergyDeposit(
    0, 0, 0, 0, posStart, posEnd, 0, 0, currEventID, 0
  );
  evNum->push_back(ID);

  // Fill RawDigit vectors from data in 'pred_fd_resps'
  auto digs = std::make_unique<std::vector<raw::RawDigit>>();

  std::vector<readout::ROPID> rIDsWritten;
  const std::string eventIDKey = "pred_fd_resps/" + std::to_string(currEventID);
  for (const std::string tpcsetStr : fFile->getGroup(eventIDKey).listObjectNames()) {
    const std::string tpcsetKey = eventIDKey + "/" + tpcsetStr;

    for (const std::string ropStr : fFile->getGroup(tpcsetKey).listObjectNames()) {
      const std::string ropKey = tpcsetKey + "/" + ropStr;

      // Get 2D response array
      std::vector<std::vector<int>> ropResp;
      HighFive::DataSet datasetRopResp = fFile->getDataSet(ropKey);
      datasetRopResp.read(ropResp);

      // Write the vector from each channel into a RawDigit vector
      const readout::ROPID rID(0, std::stoi(tpcsetStr), std::stoi(ropStr));
      for (unsigned int chNumLocal = 0; chNumLocal < ropResp.size(); chNumLocal++) {
        raw::RawDigit::ADCvector_t adcVec(fDigTickWindow);

        for (unsigned int tick = 0; tick < fDigTickWindow; tick++) {
          // Pad with zeros if desired tick window is longer
          if (tick < ropResp.size()) {
            adcVec[tick] = (short)ropResp[chNumLocal][tick];
          }
          else {
            adcVec[tick] = 0;
          }
        }

        raw::RawDigit dig(
          (raw::ChannelID_t)chNumLocal + fGeom->FirstChannelInROP(rID), adcVec.size(), adcVec
        );
        dig.SetPedestal(0);
        digs->push_back(dig);
      } // end ch

      rIDsWritten.push_back(rID);
    } // end rop
  } // end tpcset

  // Fill remaining channels with zeros
  for (const readout::ROPID rID : fGeom->Iterate<readout::ROPID>()) {
    if (std::find(rIDsWritten.begin(), rIDsWritten.end(), rID) != rIDsWritten.end()) {
      continue;
    }
    const raw::ChannelID_t firstCh = fGeom->FirstChannelInROP(rID);
    for (raw::ChannelID_t ch = firstCh; ch < firstCh + fGeom->Nchannels(rID); ch++) {
      raw::RawDigit::ADCvector_t adcVec(fDigTickWindow, 0);
      raw::RawDigit dig(ch, adcVec.size(), adcVec);
      dig.SetPedestal(0);
      digs->push_back(dig);
    }
  }

  e.put(std::move(evNum), "eventID");
  e.put(std::move(digs), "NDTranslated");
}

void extrapolation::LoadFDResp::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  std::cout << "Readng data from " << fNDFDH5FileLoc << "\n";
  fFile = new HighFive::File(fNDFDH5FileLoc, HighFive::File::ReadOnly);

  // Gather eventIDs present in file
  for (const std::string eventIDStr : fFile->getGroup("pred_fd_resps").listObjectNames()) {
    fEventIDs.push_back(std::stoi(eventIDStr));
  }
  std::sort(fEventIDs.begin(), fEventIDs.end(), std::greater<>()); // sort desc so to pop from back
  std::cout << fEventIDs.size() << " unique eventIDs in input\n";
}

void extrapolation::LoadFDResp::endJob()
{
}

DEFINE_ART_MODULE(extrapolation::LoadFDResp)
