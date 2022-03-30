#include <vector>
#include <iostream>

#include "duneextrapolation/MyCVN/func/MyCVNImageUtils.h"

cvn::MyCVNImageUtils::MyCVNImageUtils(){
  // Set a default image size
  MySetImageSize(500,500,1);
  MySetPixelMapSize(2880,500);
  // Defualt is to reverse the y view
  fViewReverse = {false,true,false};

  fUseLogScale = false;
  fDisableRegionSelection = false;
}

cvn::MyCVNImageUtils::~MyCVNImageUtils(){
}

void cvn::MyCVNImageUtils::MySetImageSize(unsigned int nWires, unsigned int nTDCs, unsigned int nViews){
  fNWires = nWires;
  fNTDCs = nTDCs;
  fNViews = nViews;
}

void cvn::MyCVNImageUtils::MySetPixelMapSize(unsigned int nWires, unsigned int nTDCs){
  fPixelMapWires = nWires;
  fPixelMapTDCs = nTDCs;
}

void cvn::MyCVNImageUtils::MySetViewReversal(bool reverseX, bool reverseY, bool reverseZ){
  fViewReverse = {reverseX,reverseY,reverseZ};
}

void cvn::MyCVNImageUtils::MySetViewReversal(std::vector<bool> reverseViews){
  if(reverseViews.size() != 3){
    std::cout << "Expected three views for view reversals... using defaults." << std::endl;
  }
  else{
    MySetViewReversal(reverseViews[0],reverseViews[1],reverseViews[2]);
  }

  return;
}

void cvn::MyCVNImageUtils::MySetLogScale(bool setLog){
  fUseLogScale = setLog;
}

void cvn::MyCVNImageUtils::ConvertPixelMapToCollectionImageVectorF(const cvn::PixelMap &pm, cvn::ImageVectorF &imageVec){
  MySetPixelMapSize(pm.fNWire,pm.fNTdc);

  // Strip out the charge vectors and use these
  std::vector<float> v0pe = pm.fPEX;
  std::vector<float> v1pe = pm.fPEY;
  std::vector<float> v2pe = pm.fPEZ;

  MyConvertChargeVectorsToImageVectorF(v0pe, v1pe, v2pe, imageVec);
}

void cvn::MyCVNImageUtils::MyConvertChargeVectorsToImageVectorF(std::vector<float> &v0pe, std::vector<float> &v1pe,
                                                           std::vector<float> &v2pe, cvn::ImageVectorF &imageVec){

  cvn::ViewVector view0;
  cvn::ViewVector view1;
  cvn::ViewVector view2;

  MyConvertChargeVectorsToViewVectors(v0pe, v1pe, v2pe, view0, view1, view2);

  // Convert the ViewVector to ViewVectorF
  cvn::ViewVectorF floatView0 = MyConvertViewVecToViewVecF(view0);
  cvn::ViewVectorF floatView1 = MyConvertViewVecToViewVecF(view1);
  cvn::ViewVectorF floatView2 = MyConvertViewVecToViewVecF(view2);

  cvn::ImageVectorF newImage = BuildCollectionImageVectorF(floatView0,floatView1,floatView2);

  imageVec = newImage;
}

void cvn::MyCVNImageUtils::MyConvertChargeVectorsToViewVectors(std::vector<float> &v0pe, std::vector<float> &v1pe, std::vector<float> &v2pe,
                                       cvn::ViewVector& view0, cvn::ViewVector& view1, cvn::ViewVector& view2){

  // Reverse requested views
  if(fViewReverse[0]) MyReverseView(v0pe);
  if(fViewReverse[1]) MyReverseView(v1pe);
  if(fViewReverse[2]) MyReverseView(v2pe);

  // Get the integrated charge for each wire
  std::vector< std::vector<float> > wireCharges;
  for (unsigned int view = 0; view < fNViews; ++view){

    std::vector<float> tempChargeVec;
    for (unsigned int wire = 0; wire < fPixelMapWires; ++wire){

      float totCharge = 0;
      for (unsigned int time = 0; time < fPixelMapTDCs; ++time){
        float val = 0.;
        unsigned int element = time + fPixelMapTDCs * wire;
        if(view == 0 ){ val = v0pe[element]; }
        if(view == 1 ){ val = v1pe[element]; }
        if(view == 2 ){ val = v2pe[element]; }
        totCharge += val;
      }
      tempChargeVec.push_back(totCharge);
    }
    wireCharges.push_back(tempChargeVec);
  }

  // Get the integrated charge for each tdc
  std::vector< std::vector<float> > tdcCharges;
  for (unsigned int view = 0; view < fNViews; ++view){

    std::vector<float> tempChargeVec;
    for (unsigned int time = 0; time < fPixelMapTDCs; ++time){


      float totCharge = 0;
      for (unsigned int wire = 0; wire < fPixelMapWires; ++wire){

        float val = 0.;
        unsigned int element = time + fPixelMapTDCs * wire;
        if(view == 0 ){ val = v0pe[element]; }
        if(view == 1 ){ val = v1pe[element]; }
        if(view == 2 ){ val = v2pe[element]; }
        totCharge += val;
      }
      tempChargeVec.push_back(totCharge);
    }
    tdcCharges.push_back(tempChargeVec);
  }

  // The output image consists of a rectangular region of the pixel map
  // We want to find the start and end wires for each view
  std::vector<unsigned int> imageStartWire(3,0);
  std::vector<unsigned int> imageEndWire(3,0);
  // And the same for TDCs
  std::vector<unsigned int> imageStartTDC(3,0);
  std::vector<unsigned int> imageEndTDC(3,0);

  if(!fDisableRegionSelection){
    // Do a rough vertex-based selection of the region for each view
    for(unsigned int view = 0; view < wireCharges.size(); ++view){
      MyGetMinMaxWires(wireCharges[view],imageStartWire[view],imageEndWire[view]);
      MyGetMinMaxTDCs(tdcCharges[view],imageStartTDC[view],imageEndTDC[view]);

//    std::cout << " Wires: " << imageStartWire[view] << ", " << imageEndWire[view] << " :: TDCs: "
//                            << imageStartTDC[view]  << ", " << imageEndTDC[view]  << std::endl;
    }
  }
  else{
    // Just use the number of wires and TDCs as the maximum values if we want to
    // use a fixed range of wires and TDC for protoDUNE's APA 3
    for(unsigned int i = 0; i < imageStartWire.size(); ++i){
      imageStartWire[i] = 0;
      imageEndWire[i] = fNWires;
      imageStartTDC[i] = 0;
      imageEndTDC[i] = fNTDCs;
    }
  }

  // Write the values to the three vectors
  for (unsigned int view = 0; view < fNViews; ++view){
    cvn::ViewVector viewChargeVec;
    for (unsigned int wire = imageStartWire[view]; wire <= imageEndWire[view]; ++wire){
      std::vector<unsigned char> wireTDCVec;
      for (unsigned int time = imageStartTDC[view]; time <= imageEndTDC[view]; ++time){

        // Get the index for the pixel map
        unsigned int element = time + fPixelMapTDCs * wire;

        // We have to convert to char and then convert back to a float
        unsigned char val = 0;
        if(view == 0){ val = MyConvertChargeToChar(v0pe[element]); }
        if(view == 1){ val = MyConvertChargeToChar(v1pe[element]); }
        if(view == 2){ val = MyConvertChargeToChar(v2pe[element]); }
        wireTDCVec.push_back(val);
      }
      viewChargeVec.push_back(wireTDCVec);
    }
    if(view == 0) view0 = viewChargeVec;
    if(view == 1) view1 = viewChargeVec;
    if(view == 2) view2 = viewChargeVec;
  }

  return;
}

cvn::ViewVectorF cvn::MyCVNImageUtils::MyConvertViewVecToViewVecF(cvn::ViewVector view){

  cvn::ViewVectorF newVec; 
  for(size_t w = 0; w < view.size(); ++w){
    std::vector<float> thisWire;
    for(size_t t = 0; t < view[w].size(); ++t){
      float chargeSC = static_cast<float>(view[w][t]);
      thisWire.push_back(chargeSC);
    }
    newVec.push_back(thisWire);
  }
  return newVec;
}

cvn::ImageVectorF cvn::MyCVNImageUtils::BuildCollectionImageVectorF(cvn::ViewVectorF v0, cvn::ViewVectorF v1, cvn::ViewVectorF v2){

  // Tensorflow wants things in the arrangement <wires, TDCs, views>
  cvn::ImageVectorF image;
  for(unsigned int w = 0; w < v2.size(); ++w){
    std::vector<std::vector<float> > wireVec;
    for(unsigned int t = 0; t < v2[0].size(); ++t){
      std::vector<float> timeVec;
      // timeVec.push_back(v0[w][t]);
      // timeVec.push_back(v1[w][t]);
      timeVec.push_back(v2[w][t]);
      wireVec.push_back(timeVec);
    } // Loop over tdcs
    image.push_back(wireVec);
  } // Loop over wires
  
  return image;
}

void cvn::MyCVNImageUtils::MyReverseView(std::vector<float> &peVec){

  std::vector<float> vecCopy(peVec.size(),0.); 

  for (unsigned int w = 0; w < fPixelMapWires; ++w)
  {
    // Get our new plane number
    unsigned int newPlane = fPixelMapWires - w - 1;

    for (unsigned int t = 0; t < fPixelMapTDCs; ++t)
    {
      float val = peVec[t + fPixelMapTDCs * w];
      vecCopy[t + fPixelMapTDCs * newPlane] = val;
    }
  }

  // Copy the values back into the original vector
  for(unsigned int e = 0; e < peVec.size(); ++e){
    float val = vecCopy[e];
    peVec[e] = val;
  }
}

void cvn::MyCVNImageUtils::MyGetMinMaxWires(std::vector<float> &wireCharges, unsigned int &minWire, unsigned int &maxWire){

  minWire = 0;
  maxWire = fNWires;

  for(unsigned int wire = 0; wire < wireCharges.size(); ++wire){

    // If we have got to fNWires from the end, the start needs to be this wire
    if(wireCharges.size() - wire == fNWires){
      break;
    }

    // For a given plane, look to see if the next 20 planes are empty. If not, this can be out start plane.
    int nEmpty = 0;
    for(unsigned int nextWire = wire + 1; nextWire <= wire + 20; ++nextWire){
      if(wireCharges[nextWire] == 0.0) ++nEmpty;
    }
    if(nEmpty < 5){
      minWire = wire;
      maxWire = wire + fNWires - 1;
      return;
    }
  }

  // If we don't find a region where we have fewer than 5 empty planes then we just want to select the fNWires wires containing
  // the most charge
  float maxCharge = 0.;
  unsigned int firstWire = 0;
  for(unsigned int wire = 0; wire < wireCharges.size() - fNWires; ++wire){
    float windowCharge = 0.;
    for(unsigned int nextwire = wire; nextwire < wire + fNWires; ++nextwire){
      windowCharge += wireCharges[nextwire];
    } 
    if(windowCharge > maxCharge){
      maxCharge = windowCharge;
      firstWire = wire;
    }
  } 
  minWire = firstWire;
  maxWire = firstWire + fNWires - 1;

  std::cout << "Used alternate method to get min and max wires due to vertex determination failure: " << minWire << ", " << maxWire << std::endl;
}

void cvn::MyCVNImageUtils::MyGetMinMaxTDCs(std::vector<float> &tdcCharges, unsigned int &minTDC, unsigned int &maxTDC){

  minTDC = 0;
  maxTDC = fNTDCs;

  for(unsigned int tdc = 0; tdc < tdcCharges.size(); ++tdc){

    // If we have got to fNWires from the end, the start needs to be this wire
    if(tdcCharges.size() - tdc == fNTDCs){
      break;
    } 

    // For a given tdc, look to see if the next 20 tdcs are empty. If not, this can be out start tdc.
    int nEmpty = 0;
    for(unsigned int nextTDC = tdc + 1; nextTDC <= tdc + 20; ++nextTDC){
      if(tdcCharges[nextTDC] == 0.0) ++nEmpty;
    } 
    if(nEmpty < 5){
      minTDC = tdc;
      maxTDC = tdc + fNTDCs - 1;
      return;
    } 
  }

  // If we don't find a region where we have fewer than 5 empty tdcs then we just want to select the fNTDCs tdcs containing
  // the most charge
  float maxCharge = 0.;
  unsigned int firstTDC = 0;
  for(unsigned int tdc = 0; tdc < tdcCharges.size() - fNTDCs; ++tdc){
    float windowCharge = 0.;
    for(unsigned int nexttdc = tdc; nexttdc < tdc + fNTDCs; ++nexttdc){
      windowCharge += tdcCharges[nexttdc];
    }
    if(windowCharge > maxCharge){
      maxCharge = windowCharge;
      firstTDC = tdc;
    }
  }
  minTDC = firstTDC;
  maxTDC = firstTDC + fNTDCs - 1;

  std::cout << "Used alternate method to get min and max tdcs due to vertex determination failure: " << minTDC << ", " << maxTDC << std::endl;
}

unsigned char cvn::MyCVNImageUtils::MyConvertChargeToChar(float charge){

  float peCorrChunk;
  float truncateCorr;
  float centreScale = 0.7;
  if(fUseLogScale){
    float scaleFrac=(log(charge)/log(1000));
    truncateCorr= ceil(centreScale*scaleFrac*255.0);
      }
  else{
    peCorrChunk = (1000.) / 255.0;
    truncateCorr = ceil((charge)/(peCorrChunk));
  }
  if (truncateCorr > 255) return (unsigned char)255;
  else return (unsigned char)truncateCorr;
}