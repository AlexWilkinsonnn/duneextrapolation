////////////////////////////////////////////////////////////////////////
// Class:       AddFDResp
// Plugin Type: analyzer (Unknown Unknown)
// File:        AddFDResp_module.cc
//
// Crated on 22 Nov 23 Alex Wilkinson
// Reads in 3d packets from ND-FD pair HDF5 file. Translates and ECC
// rotates ND event s.t. it is align with FD event. Project ND packets
// to (channel, tick). Write out to ND-FD pair HDF5 the projected ND
// packets and FD wire response for relevant TPCs
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
#include "dunereco/CVN/func/Result.h"
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
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
#include <map>
#include <math.h>

typedef struct packet3d {
  int eventID;
  int adc;
  double x;
  double x_module;
  double y;
  double z;
  double z_module;
  unsigned int forward_facing_anode;
} packet3d;

HighFive::CompoundType make_packet3d() {
  return {
    {"eventID", HighFive::AtomicType<int>{}},
    {"adc", HighFive::AtomicType<int>{}},
    {"x", HighFive::AtomicType<double>{}},
    {"x_module", HighFive::AtomicType<double>{}},
    {"y", HighFive::AtomicType<double>{}},
    {"z", HighFive::AtomicType<double>{}},
    {"z_module", HighFive::AtomicType<double>{}},
    {"forward_facing_anode", HighFive::AtomicType<unsigned int>{}}
  };
}

typedef struct vertex {
  int eventID;
  double x_vert;
  double y_vert;
  double z_vert;
} vertex;

HighFive::CompoundType make_vertex() {
  return {
    {"eventID", HighFive::AtomicType<int>{}},
    {"x_vert", HighFive::AtomicType<double>{}},
    {"y_vert", HighFive::AtomicType<double>{}},
    {"z_vert", HighFive::AtomicType<double>{}}
  };
}

typedef struct packetProj {
  int adc;
  int local_ch;
  int tick;
  double nd_drift_dist;
  double fd_drift_dist;
  double nd_x_module;
  double wire_dist;
  int forward_facing_anode;
} packetProj;

HighFive::CompoundType make_packetProj() {
  return {
    {"adc", HighFive::AtomicType<int>{}},
    {"local_ch", HighFive::AtomicType<int>{}},
    {"tick", HighFive::AtomicType<int>{}},
    {"nd_drift_dist", HighFive::AtomicType<double>{}},
    {"fd_drift_dist", HighFive::AtomicType<double>{}},
    {"nd_x_module", HighFive::AtomicType<double>{}},
    {"wire_dist", HighFive::AtomicType<double>{}},
    {"forward_facing_anode", HighFive::AtomicType<int>{}}
  };
}

typedef struct vertexProj {
  int eventID;
  int local_ch;
  int tick;
  int tpc_set;
  int readout;
} vertexProj;

HighFive::CompoundType make_vertexProj() {
  return {
    {"eventID", HighFive::AtomicType<int>{}},
    {"local_ch", HighFive::AtomicType<int>{}},
    {"tick", HighFive::AtomicType<int>{}},
    {"tpc_set", HighFive::AtomicType<int>{}},
    {"readout", HighFive::AtomicType<int>{}}
  };
}

HIGHFIVE_REGISTER_TYPE(packet3d, make_packet3d)
HIGHFIVE_REGISTER_TYPE(vertex, make_vertex)
HIGHFIVE_REGISTER_TYPE(packetProj, make_packetProj)
HIGHFIVE_REGISTER_TYPE(vertexProj, make_vertexProj)

namespace extrapolation {
  class AddFDResp;
}

class extrapolation::AddFDResp : public art::EDAnalyzer {
public:
  explicit AddFDResp(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  AddFDResp(AddFDResp const&) = delete;
  AddFDResp(AddFDResp&&) = delete;
  AddFDResp& operator=(AddFDResp const&) = delete;
  AddFDResp& operator=(AddFDResp&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Methods
  bool inWireCellBoundingBox(const double x, const double y, const double z);
  void alignNDWithFD(std::vector<packet3d> &packets, const vertex& NDVtx, const vertex& FDVtx);

  // Members
  const geo::GeometryCore* fGeom;

  HighFive::File* fFile;

  // Data to write out to hdf5
  std::vector<packetProj> fPacketProjsZ;
  std::vector<packetProj> fPacketProjsU;
  std::vector<packetProj> fPacketProjsV;

  // Data to read in from hdf5
  std::map<int, vertex> fNDVertices;
  std::map<int, vertex> fFDVertices;
  std::map<int, std::vector<packet3d>> fNDPackets;

  // Compound data to write out to hdf5
  std::vector<vertexProj> fNDFDVerticesProj;

  // Product labels
  std::string fEventIDSEDLabel;
  std::string fFDSEDLabel;
  std::string fRawDigitLabel;

  std::string fNDFDH5FileLoc;
  double fECCRotation;
  std::vector<std::vector<std::vector<double>>> fWireCellAPABoundingBoxes;
  double fNDProjForwardAnodeXShift;
  double fNDProjBackwardAnodeXShift;
  bool fProjectFDDepos;
  unsigned int fMaxTick;
};


extrapolation::AddFDResp::AddFDResp(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fEventIDSEDLabel           (p.get<std::string>("EventIDSEDLabel")),
    fFDSEDLabel                (p.get<std::string>("FDSEDLabel")),
    fRawDigitLabel             (p.get<std::string>("RawDigitLabel")),
    fNDFDH5FileLoc             (p.get<std::string>("NDFDH5FileLoc")),
    fECCRotation               (p.get<double>("ECCRotation")),
    fWireCellAPABoundingBoxes  (p.get<std::vector<std::vector<std::vector<double>>>>("WireCellAPABoundingBoxes")),
    fNDProjForwardAnodeXShift  (p.get<double>("NDProjForwardAnodeXShift")),
    fNDProjBackwardAnodeXShift (p.get<double>("NDProjBackwardAnodeXShift")),
    fProjectFDDepos            (p.get<bool>("ProjectNDDepos")),
    fMaxTick                   (p.get<unsigned int>("MaxTick"))
{
  consumes<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);
  consumes<std::vector<raw::RawDigit>>(fRawDigitLabel);

  if (fProjectFDDepos) {
    consumes<std::vector<sim::SimEnergyDeposit>>(fFDSEDLabel);
  }

  fMaxTick = std::min(fMaxTick, (unsigned int)6000);
}

void extrapolation::AddFDResp::analyze(art::Event const& e)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  // Get eventID
  const auto eventIDSED = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);
  int eventID = (*eventIDSED)[0].TrackID();

  // Get RawDigits
  const auto digits = e.getValidHandle<std::vector<raw::RawDigit>>(fRawDigitLabel);

  // Get vertices and depos from read from HDF5
  const vertex NDVtx = fNDVertices[eventID];
  const vertex FDVtx = fFDVertices[eventID];
  std::vector<packet3d> packets = fNDPackets[eventID];

  alignNDWithFD(packets, NDVtx, FDVtx);

  // Get packet wire projections
  std::map<readout::ROPID, std::vector<packetProj>> eventPacketProjs;
  for (const auto packet : packets) {
    // Discard packets not in wirecell bounding boxes for drifitng
    if (!inWireCellBoundingBox(packet.x, packet.y, packet.z)) {
      continue;
    }

    const geo::Point_t packetLoc(packet.x, packet.y, packet.z);

    // Make projections to all relevant readouts (packet would be drifted to two induction and the
    // collection plane that faces it)
    const geo::TPCID tID = fGeom->PositionToTPCID(packetLoc);
    for (const geo::PlaneID pID : fGeom->Iterate<geo::PlaneID>(tID)) {
      const geo::PlaneGeo pGeo = fGeom->Plane(pID);
      const readout::ROPID rID = fGeom->WirePlaneToROP(pID);
    
      const raw::ChannelID_t ch = 
        fGeom->NearestChannel(packetLoc, pID) - fGeom->FirstChannelInROP(rID);

      double xShift = 0.0;
      if (packet.forward_facing_anode == 1) {
        xShift += fNDProjForwardAnodeXShift;
      }
      else if (packet.forward_facing_anode == 0) {
        xShift += fNDProjBackwardAnodeXShift;
      }
      else {
        throw cet::exception("AddFDResp")
          << "Packet has unset forward_facing_anode"
          << " - Line " << __LINE__ << " in file " << __FILE__ << "\n";
      }

      const int tick = (int)detProp.ConvertXToTicks(packet.x + xShift, pID);
      if (tick >= (int)fMaxTick) {
        continue;
      }

      const float driftDistanceFD = pGeo.DistanceFromPlane(packetLoc);

      // Gets the distance in the direction perpendicular to wire to the closest wire
      const double wireCoord = pGeo.WireCoordinate(packetLoc);
      const double wireDistance = (wireCoord - (double)(int)(0.5 + wireCoord)) * pGeo.WirePitch();

      const packetProj p = {
        packet.adc,
        (int)ch,
        tick,
        packet.z_module,
        driftDistanceFD,
        packet.x_module,
        wireDistance,
        (int)packet.forward_facing_anode
      };
      eventPacketProjs[rID].push_back(p);
    }
  }

  // Write projections to the HDF5
  for (const auto rID_projs : eventPacketProjs) {
    const readout::ROPID rID = rID_projs.first;
    const std::vector<packetProj> projs = rID_projs.second;

    const std::string groupPath =
      "nd_packet_wire_projs/" +
      std::to_string(eventID) + "/" +
      std::to_string(rID.TPCset) + "/" +
      std::to_string(rID.ROP);
    fFile->createDataSet(groupPath, projs);
  }

  // Find the projection for the ND-FD aligned vertex (if there is a valid one)
  if (!inWireCellBoundingBox(FDVtx.x_vert, FDVtx.y_vert, FDVtx.z_vert)) {
    const vertexProj FDVtxProj = { eventID, -1, -1, -1, -1 };
    fNDFDVerticesProj.push_back(FDVtxProj);
  }
  else {
    const geo::Point_t vtxLoc(FDVtx.x_vert, FDVtx.y_vert, FDVtx.z_vert);

    const geo::TPCID tID = fGeom->PositionToTPCID(vtxLoc);
    for (const geo::PlaneID pID : fGeom->Iterate<geo::PlaneID>(tID)) {
      const readout::ROPID rID = fGeom->WirePlaneToROP(pID);

      const raw::ChannelID_t ch = 
        fGeom->NearestChannel(vtxLoc, pID) - fGeom->FirstChannelInROP(rID);

      const int tick = (int)detProp.ConvertXToTicks(FDVtx.x_vert, pID);

      const vertexProj FDVtxProj = { eventID, (int)ch, tick, (int)rID.TPCset, (int)rID.ROP };
      fNDFDVerticesProj.push_back(FDVtxProj);
    }
  }

  // Read in relevant RawDigits
  std::map<readout::ROPID, std::vector<std::vector<short>>> eventRawDigits;
  for (const raw::RawDigit& dig : *digits) { 
    const readout::ROPID rID = fGeom->ChannelToROP(dig.Channel());

    // Skip if no ND packets for this ROP
    if (eventPacketProjs.find(rID) == eventPacketProjs.end()) {
      continue;
    }

    if (eventRawDigits.find(rID) == eventRawDigits.end()) {
      eventRawDigits[rID] = 
        std::vector<std::vector<short>>(fGeom->Nchannels(rID), std::vector<short>(fMaxTick, 0));
    }

    raw::RawDigit::ADCvector_t adcs(dig.Samples());
    raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

    for (unsigned int tick = 0; tick < fMaxTick; tick++) { 
      const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;
      eventRawDigits[rID][dig.Channel() - fGeom->FirstChannelInROP(rID)][tick] = adc;
    }
  }

  // Write RawDigits to the HDF5
  for (const auto rID_adcs : eventRawDigits) { 
    const readout::ROPID rID = rID_adcs.first;
    const std::vector<std::vector<short>> adcs = rID_adcs.second;

    const std::string groupPath =
      "fd_resp/" +
      std::to_string(eventID) + "/" +
      std::to_string(rID.TPCset) + "/" +
      std::to_string(rID.ROP);
    // Want compression for these
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{100, 750}));
    props.add(HighFive::Deflate(9));
    fFile->createDataSet(groupPath, adcs, props);
  }

  if (fProjectFDDepos) {
    std::map<readout::ROPID, std::vector<packetProj>> eventFDSEDProjs;
    const auto SEDs = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fFDSEDLabel);

    for (const auto SED : *SEDs) { 
      if (!inWireCellBoundingBox(SED.X(), SED.Y(), SED.Z())) {
        continue;
      }

      const geo::Point_t SEDLoc = SED.MidPoint();

      const geo::TPCID tID = fGeom->PositionToTPCID(SEDLoc);
      for (const geo::PlaneID pID : fGeom->Iterate<geo::PlaneID>(tID)) {
        const geo::PlaneGeo pGeo = fGeom->Plane(pID);
        const readout::ROPID rID = fGeom->WirePlaneToROP(pID);
      
        const raw::ChannelID_t ch = 
          fGeom->NearestChannel(SEDLoc, pID) - fGeom->FirstChannelInROP(rID);

        const int tick = (int)detProp.ConvertXToTicks(SED.X(), pID);

        const float driftDistanceFD = pGeo.DistanceFromPlane(SEDLoc);

        // Gets the distance in the direction perpendicular to wire to the closest wire
        const double wireCoord = pGeo.WireCoordinate(SEDLoc);
        const double wireDistance = (wireCoord - (double)(int)(0.5 + wireCoord)) * pGeo.WirePitch();

        const packetProj p = {
          SED.NumElectrons(),
          (int)ch,
          tick,
          -1,
          driftDistanceFD,
          -1,
          wireDistance,
          -1
        };
        eventFDSEDProjs[rID].push_back(p);
      }
    }

    for (const auto rID_projs : eventFDSEDProjs) {
      const readout::ROPID rID = rID_projs.first;
      const std::vector<packetProj> projs = rID_projs.second;

      if (eventPacketProjs.find(rID) == eventPacketProjs.end()) {
        continue;
      }

      const std::string groupPath =
        "fd_deps_wire_projs/" +
        std::to_string(eventID) + "/" +
        std::to_string(rID.TPCset) + "/" +
        std::to_string(rID.ROP);
      fFile->createDataSet(groupPath, projs);
    }
  }
}

void extrapolation::AddFDResp::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  fFile = new HighFive::File(fNDFDH5FileLoc, HighFive::File::ReadWrite);

  // Read in depos and vertices
  std::vector<vertex> FDVtxs;
  HighFive::DataSet datasetFDVtxs = fFile->getDataSet("fd_vertices");
  datasetFDVtxs.read(FDVtxs);
  std::vector<vertex> NDVtxs;
  HighFive::DataSet datasetNDVtxs = fFile->getDataSet("vertices");
  datasetNDVtxs.read(NDVtxs);
  std::vector<packet3d> NDPackets;
  HighFive::DataSet datasetNDPackets = fFile->getDataSet("3d_packets");
  datasetNDPackets.read(NDPackets);

  for (std::size_t i = 0; i < NDPackets.size(); i++) {
    // ND packets use z for drift coordinate, swap this here to align with FD things
    const double x = NDPackets[i].x;
    const double z = NDPackets[i].z;
    NDPackets[i].x = z;
    NDPackets[i].z = x;
    fNDPackets[NDPackets[i].eventID].push_back(NDPackets[i]);
  }
  for (const vertex NDVtx : NDVtxs) {
    fNDVertices[NDVtx.eventID] = NDVtx;
  }
  for (const vertex FDVtx : FDVtxs) {
    fFDVertices[FDVtx.eventID] = FDVtx;
  }
}

void extrapolation::AddFDResp::endJob()
{
  // Write out the projected ND vertices
  fFile->createDataSet("ndfd_vertices_projs", fNDFDVerticesProj);
}

void extrapolation::AddFDResp::alignNDWithFD(
  std::vector<packet3d> &packets, const vertex& NDVtx, const vertex& FDVtx
)
{
  // Shift the packets s.t. vertex is (0,0,0) to prepare for ECC rotation
  for (std::size_t i = 0; i < packets.size(); i++) {
    packets[i].x = packets[i].x - NDVtx.x_vert;
    packets[i].y = packets[i].y - NDVtx.y_vert;
    packets[i].z = packets[i].z - NDVtx.z_vert;
  }

  // Apply ECC rotation to packets to fully align ND and FD responses
  // ECC rotation should be -0.202 rad (clockwise) about x (drift direction)
  for (std::size_t i = 0; i < packets.size(); i++) {
    const double y = packets[i].y;
    const double z = packets[i].z;
    packets[i].y = y * cos(fECCRotation) - z * sin(fECCRotation);
    packets[i].z = y * sin(fECCRotation) + z * cos(fECCRotation);
  }

  // Shift the packets s.t. the ND and FD vertices are align
  for (std::size_t i = 0; i < packets.size(); i++) {
    packets[i].x = packets[i].x + FDVtx.x_vert;
    packets[i].y = packets[i].y + FDVtx.y_vert;
    packets[i].z = packets[i].z + FDVtx.z_vert;
  }
}

bool extrapolation::AddFDResp::inWireCellBoundingBox(
  const double x, const double y, const double z
)
{
  for (const std::vector<std::vector<double>> range : fWireCellAPABoundingBoxes) {
    if (
      x >= range[0][0] && x <= range[1][0] &&
      y >= range[0][1] && y <= range[1][1] &&
      z >= range[0][2] && z <= range[1][2]
    ) {
      return true;
    }
  }
  return false;
}

DEFINE_ART_MODULE(extrapolation::AddFDResp)
