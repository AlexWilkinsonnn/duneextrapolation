////////////////////////////////////////////////////////////////////////
// Class:       AddNDProj
// Plugin Type: analyzer (Unknown Unknown)
// File:        AddNDProj_module.cc
//
// Crated on 06 Nov 24 Alex Wilkinson
// Reads in 3d packets infilled from ND-FD pair HDF5 file.
// Translates and ECC rotates ND event s.t. it is align with FD event.
// Project ND packets to (channel, tick).
// The file from LoadFDDepos is still used for the eventID and
// getting the detector properties because I am dumb.
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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

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
  int adc; // NOTE better as float, used int in training dataset so keep as int for now
  double x;
  double x_module;
  double y;
  double z;
  double z_module;
  int forward_facing_anode;
  int infilled;
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
    {"forward_facing_anode", HighFive::AtomicType<int>{}},
    {"infilled", HighFive::AtomicType<int>{}}
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
  int infilled;
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
    {"forward_facing_anode", HighFive::AtomicType<int>{}},
    {"infilled", HighFive::AtomicType<int>{}}
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
  class AddNDProj;
}

class extrapolation::AddNDProj : public art::EDAnalyzer {
public:
  explicit AddNDProj(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  AddNDProj(AddNDProj const&) = delete;
  AddNDProj(AddNDProj&&) = delete;
  AddNDProj& operator=(AddNDProj const&) = delete;
  AddNDProj& operator=(AddNDProj&&) = delete;

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

  // See config fcl for description of these members
  std::string fNDFDH5FileLoc;
  double fECCRotation;
  std::vector<std::vector<std::vector<double>>> fWireCellAPABoundingBoxes;
  double fNDProjForwardAnodeXShift;
  double fNDProjBackwardAnodeXShift;
  unsigned int fMaxTick;
  std::string fPacketsH5DataSetName;
};


extrapolation::AddNDProj::AddNDProj(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fEventIDSEDLabel           (p.get<std::string>("EventIDSEDLabel")),
    fNDFDH5FileLoc             (p.get<std::string>("NDFDH5FileLoc")),
    fECCRotation               (p.get<double>("ECCRotation")),
    fWireCellAPABoundingBoxes  (p.get<std::vector<std::vector<std::vector<double>>>>("WireCellAPABoundingBoxes")),
    fNDProjForwardAnodeXShift  (p.get<double>("NDProjForwardAnodeXShift")),
    fNDProjBackwardAnodeXShift (p.get<double>("NDProjBackwardAnodeXShift")),
    fMaxTick                   (p.get<unsigned int>("MaxTick")),
    fPacketsH5DataSetName      (p.get<std::string>("PacketsH5DataSetName"))
{
  consumes<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);

  fMaxTick = std::min(fMaxTick, (unsigned int)6000);
}

void extrapolation::AddNDProj::analyze(art::Event const& e)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  // Get eventID
  const auto eventIDSED = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEventIDSEDLabel);
  int eventID = (*eventIDSED)[0].TrackID();

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
      else if (packet.infilled == 0) { // Only expected to be valid for non-infilled
        throw cet::exception("AddNDProj")
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
        packet.forward_facing_anode,
        packet.infilled
      };
      eventPacketProjs[rID].push_back(p);
    }
  }

  // Write projections to the HDF5
  for (const auto& rID_projs : eventPacketProjs) {
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
}

void extrapolation::AddNDProj::beginJob()
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
  HighFive::DataSet datasetNDPackets = fFile->getDataSet(fPacketsH5DataSetName);
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

void extrapolation::AddNDProj::endJob()
{
  // Write out the projected ND vertices
  fFile->createDataSet("ndfd_vertices_projs", fNDFDVerticesProj);
}

void extrapolation::AddNDProj::alignNDWithFD(
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

bool extrapolation::AddNDProj::inWireCellBoundingBox(
  const double x, const double y, const double z
)
{
  for (const std::vector<std::vector<double>>& range : fWireCellAPABoundingBoxes) {
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

DEFINE_ART_MODULE(extrapolation::AddNDProj)
