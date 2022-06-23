////////////////////////////////////////////////////////////////////////
// Class: ProjectionData, PacketData
// File:  Types.h
//
// Created on Thu 23 Jun 22 by Alex Wilkinson.
////////////////////////////////////////////////////////////////////////

#ifndef EXTRAPOLATION_TYPES_H
#define EXTRAPOLATION_TYPES_H

#include "larcorealg/Geometry/GeometryCore.h"

namespace extrapolation
{

  struct PacketData
  {
    PacketData(geo::Point_t _pos, int _adc, double _NDDrift)
      : pos(_pos), adc(_adc), NDDrift(_NDDrift)
    { }

    geo::Point_t pos;
    int          adc;
    double       NDDrift;
  };

  struct ProjectionData
  {
    ProjectionData(PacketData& _packet, int _localCh, int _tick, double _FDDrift,
       double _wireDistance)
      : packet(_packet), localCh(_localCh), tick(_tick), FDDrift(_FDDrift),
        wireDistance(_wireDistance)
    { }
    ProjectionData(PacketData& _packet)
      : packet(_packet)
    { }

    PacketData&  packet;
    int          localCh;
    int          tick;
    int          FDDrift;
    double       wireDistance;
  };

}

#endif
