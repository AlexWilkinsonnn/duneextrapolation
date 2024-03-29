////////////////////////////////////////////////////////////////////////
/// \file    TFNetHandler.cxx
/// \brief   TFNetHandler for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
///          Leigh Whitehead   - leigh.howard.whitehead@cern.ch
///          Saul Alonso Monsalve - saul.alonso.monsalve@cern.ch
///
/// Copied from dunereco/CVN/art on 15 Mar 2022
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <string>
#include "cetlib/getenv.h"

#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "duneextrapolation/Misc/MyCVN/art/MyTFNetHandler.h"
#include "duneextrapolation/Misc/MyCVN/func/MyCVNImageUtils.h"

namespace cvn
{

  MyTFNetHandler::MyTFNetHandler(const fhicl::ParameterSet& pset):
    fLibPath(cet::getenv(pset.get<std::string>("LibPath", ""))),
    fTFProtoBuf  (fLibPath+"/"+pset.get<std::string>("TFProtoBuf")),
    fUseLogChargeScale(pset.get<bool>("ChargeLogScale")),
    fImageWires(pset.get<unsigned int>("NImageWires")),
    fImageTDCs(pset.get<unsigned int>("NImageTDCs")),
    fReverseViews(pset.get<std::vector<bool> >("ReverseViews"))
  {

    // Construct the TF Graph object. The empty vector {} is used since the protobuf
    // file gives the names of the output layer nodes
    mf::LogInfo("MyTFNetHandler") << "Loading network: " << fTFProtoBuf << std::endl;
    // const std::vector<std::string> outputs = {"Identity_0", "Identity_1", "Identity_2", 
    //   "Identity_3", "Identity_4", "Identity_5","Identity_6"};
    const std::vector<std::string> outputs = {"resnext/is_antineutrino/Sigmoid", "resnext/flavour/Softmax",
      "resnext/interaction/Softmax", "resnext/protons/Exp", "resnext/pions/Exp", "resnext/pizeros/Exp",
      "resnext/neutrons/Exp"};
    fTFGraph = tf::MyGraph::Mycreate(fTFProtoBuf.c_str(),outputs,pset.get<int>("NInputs"),pset.get<int>("NOutputs"));
    if(!fTFGraph){
      art::Exception(art::errors::Unknown) << "Tensorflow model not found or incorrect";
    }

  }
 
  // Check the network outputs
  bool check(const std::vector< std::vector< float > > & outputs)
  {
    if (outputs.size() == 1) return true;
    size_t aux = 0;
    for (size_t o = 0; o < outputs.size(); ++o)
    {   
        size_t aux2 = 0;

        for (size_t i = 0; i < outputs[o].size(); ++i)
            if (outputs[o][i] == 0.0 || outputs[o][i] == 1.0)
                aux2++;
        if (aux2 == outputs[o].size()) aux++;
    }
    return aux == outputs.size() ? false : true;        
  }

  // Fill outputs with value -3
  void fillEmpty(std::vector< std::vector< float > > & outputs)
  {
    for (size_t o = 0; o < outputs.size(); ++o)
    {
        for (size_t i = 0; i < outputs[o].size(); ++i)
            outputs[o][i] = -3.0;
    }
    return;
  }

  std::vector< std::vector<float> > MyTFNetHandler::Predict(const PixelMap& pm)
  {
   
    MyCVNImageUtils imageUtils;

    // Configure the image utility  
    imageUtils.MySetViewReversal(fReverseViews);
    imageUtils.MySetImageSize(fImageWires,fImageTDCs,3);
    imageUtils.MySetLogScale(fUseLogChargeScale);

    ImageVectorF thisImage;
    // imageUtils.ConvertPixelMapToImageVectorF(pm,thisImage);
    imageUtils.ConvertPixelMapToCollectionImageVectorF(pm, thisImage);
    std::vector<ImageVectorF> vecForTF;

    // for (auto vec1 : thisImage) {
    //   for (auto vec2 : vec1) { 
    //     for (auto val : vec2) {
    //       if (val != 0) {
    //         std::cout << val << " ";
    //       }
    //     }
    //   }
    // }
    // std::cout << "\n";

    vecForTF.push_back(thisImage);

    std::vector< std::vector< std::vector< float > > > cvnResults; // shape(samples, #outputs, output_size)
    bool status = false;

    int counter = 0;

    do{ // do until it gets a correct result
        // std::cout << "Number of CVN result vectors " << cvnResults.size() << " with " << cvnResults[0].size() << " categories" << std::endl;
        cvnResults = fTFGraph->Myrun(vecForTF);
        status = check(cvnResults[0]);
        //std::cout << "Status: " << status << std::endl;
        counter++;
        if(counter==10){
            std::cout << "Error, CVN never outputing a correct result. Filling result with zeros.";
            std::cout << std::endl;
            fillEmpty(cvnResults[0]);
            break;
        }
    }while(status == false);

    // Convert proton/pion/pizero/neutron number regression outputs to the expected CVN class outputs.
    for (int i = 3; i < 7; i++) {
      if ((int)(cvnResults[0][i][0] + 0.5) == 0) {
        cvnResults[0][i] = {1.0, 0.0, 0.0, 0.0};
      }
      else if ((int)(cvnResults[0][i][0] + 0.5) == 1) {
        cvnResults[0][i] = {0.0, 1.0, 0.0, 0.0};
      }
      else if ((int)(cvnResults[0][i][0] + 0.5) == 2) {
        cvnResults[0][i] = {0.0, 0.0, 1.0, 0.0};
      }
      else if ((int)(cvnResults[0][i][0] + 0.5) > 2) {
        cvnResults[0][i] = {0.0, 0.0, 0.0, 1.0};
      }
    }

    std::cout << "Classifier summary: ";
    std::cout << std::endl;
    int output_index = 0;
    for(auto const & output : cvnResults[0])
    {
      std::cout << "Output " << output_index++ << ": ";
      for(auto const v : output)
          std::cout << v << ", ";
      std::cout << std::endl;
    }
    std::cout << std::endl;

    return cvnResults[0];
  }

  /* 
  // The standard output has 13 elements, this function sums the convenient ones 
  std::vector<float> MyTFNetHandler::PredictFlavour(const PixelMap& pm){

    std::vector<float> fullResults = this->Predict(pm);

    std::vector<float> flavourResults;

    // First element is CC numu
    float sumNumu  = fullResults[0] + fullResults[1] + fullResults[2] + fullResults[3];
    // Then CC nue
    float sumNue   = fullResults[4] + fullResults[5] + fullResults[6] + fullResults[7];
    // Then CC nutau
    float sumNutau = fullResults[8] + fullResults[9] + fullResults[10] + fullResults[11];
    // End with NC
    float sumNC    = fullResults[12];

    flavourResults.push_back(sumNumu);
    flavourResults.push_back(sumNue);
    flavourResults.push_back(sumNutau);
    flavourResults.push_back(sumNC);

    return flavourResults;
  }
  */

}

