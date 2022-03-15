////////////////////////////////////////////////////////////////////////
/// \file    CVNImageUtils.h
/// \brief   Utilities for producing images for the CVN
/// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
///
/// Copied from dunereco/CVN/func on 15 Mar 2022.
////////////////////////////////////////////////////////////////////////

#ifndef MYCVN_IMAGE_UTILS_H
#define MYCVN_IMAGE_UTILS_H

#include <vector>

#include "dunereco/CVN/func/PixelMap.h"

namespace cvn
{

  /// Useful typedefs
  typedef std::vector<std::vector<unsigned char> > ViewVector;
  typedef std::vector<ViewVector> ImageVector;
  typedef std::vector<std::vector<float> > ViewVectorF;
  typedef std::vector<ViewVectorF> ImageVectorF;

  /// Class containing some utility functions for all things CVN
  class MyCVNImageUtils
  {
  public:
    MyCVNImageUtils();
    ~MyCVNImageUtils();

    /// Function to set any views that need reversing
    void MySetViewReversal(bool reverseX, bool reverseY, bool reverseZ);
    void MySetViewReversal(std::vector<bool> reverseViews);

    /// Set the log scale for charge
    void MySetLogScale(bool setLog);

    /// Set up the image size that we want to have
    void MySetImageSize(unsigned int nWires, unsigned int nTDCs, unsigned int nViews);

    /// Set the input pixel map size
    void MySetPixelMapSize(unsigned int nWires, unsigned int nTDCs);

    /// Convert the hit charge into the range 0 to 255 required by the CVN
    unsigned char MyConvertChargeToChar(float charge);

    /// Convert a pixel map into an image vector (float version)
    void ConvertPixelMapToCollectionImageVectorF(const PixelMap &pm, ImageVectorF &imageVec);

    /// Float version of conversion for convenience of TF interface
    void MyConvertChargeVectorsToImageVectorF(std::vector<float> &v0pe, std::vector<float> &v1pe,
                                           std::vector<float> &v2pe, ImageVectorF &imageVec);  

  private:

    /// Convert a ViewVector into a ViewVectorF
    ViewVectorF MyConvertViewVecToViewVecF(ViewVector view);

    /// Base function for conversion of the Pixel Map to our required output format
    void MyConvertChargeVectorsToViewVectors(std::vector<float> &v0pe, std::vector<float> &v1pe, std::vector<float> &v2pe,
                                  ViewVector& view0, ViewVector& view1, ViewVector& view2);

    /// Make the image vector from the view vectors
    ImageVectorF BuildCollectionImageVectorF(ViewVectorF v0, ViewVectorF v1, ViewVectorF v2);

    /// Get the minimum and maximum wires from the pixel map needed to make the image
    void MyGetMinMaxWires(std::vector<float> &wireCharges, unsigned int &minWire, unsigned int &maxWire); 

    /// Get the minimum and maximum tdcs from the pixel map needed to make the image
    void MyGetMinMaxTDCs(std::vector<float> &tdcCharges, unsigned int &minTDC, unsigned int &maxTDC); 

    /// Funtion to actually reverse the view
    void MyReverseView(std::vector<float> &peVec);

    /// Number of views of each event
    unsigned int fNViews;

    /// Number of wires to use for the image width
    unsigned int fNWires;

    /// Number of TDCs to use for the image height
    unsigned int fNTDCs;

    /// Input pixel map sizes
    unsigned int fPixelMapWires;
    unsigned int fPixelMapTDCs;

    /// Vector of bools to decide if any views need to be reversed
    std::vector<bool> fViewReverse;

    /// Disable the region finding?
    bool fDisableRegionSelection;

    /// Use a log scale for charge?
    bool fUseLogScale;

  };

}

#endif  // MYCVN_IMAGE_UTILS_H
