////////////////////////////////////////////////////////////////////////
/// \file    TFNetHandler.h
/// \brief   TFNetHandler for CVN
/// \author  Leigh Whitehead
///
/// Copied from dunereco/CVN/art on 15 Mar 2022.
////////////////////////////////////////////////////////////////////////

#ifndef CVN_MYTFNETHANDLER_H
#define CVN_MYTFNETHANDLER_H

#include <vector>
#include <memory>

#include "dunereco/CVN/func/PixelMap.h"
#include "dunereco/CVN/func/InteractionType.h"
#include "fhiclcpp/ParameterSet.h"
#include "duneextrapolation/Misc/MyCVN/tf/Mytf_graph.h"

namespace cvn
{

  /// Wrapper for caffe::Net which handles construction and prediction
  class MyTFNetHandler
  {
  public:

    /// Constructor which takes a pset with DeployProto and ModelFile fields
    MyTFNetHandler(const fhicl::ParameterSet& pset);

    /// Number of outputs in neural net
    int NOutput() const;

    /// Number of outputs in neural net
    int NFeatures() const;

    /// Return prediction arrays for PixelMap
    std::vector< std::vector<float> > Predict(const PixelMap& pm);

    /// Return four element vector with summed numu, nue, nutau and NC elements
    std::vector<float> PredictFlavour(const PixelMap& pm);

  private:

    std::string  fLibPath;  ///< Library path (typically dune_pardata...)
    std::string  fTFProtoBuf;  ///< location of the tf .pb file in the above path
    bool         fUseLogChargeScale;  ///< Is the charge using a log scale?
    unsigned int fImageWires;  ///< Number of wires for the network to classify
    unsigned int fImageTDCs;   ///< Number of tdcs for the network to classify
    std::vector<bool> fReverseViews; ///< Do we need to reverse any views?
    std::unique_ptr<tf::MyGraph> fTFGraph; ///< Tensorflow graph

  };

}

#endif  // CVN_MyTFNetHandler_H
