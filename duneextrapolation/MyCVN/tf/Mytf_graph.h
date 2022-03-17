////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class:       Graph
//// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
///               P.Plonski,                      from DUNE, WUT, Sept. 2017
////              S. Alonso Monsalve,             from DUNE, CERN, Aug. 2018
//// Iterface to run Tensorflow graph saved to a file. First attempts, almost functional.
////
/// Copied from dunereco/CVN/tf on 16 Mar 2022
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MyGraph_h
#define MyGraph_h

#include <memory>
#include <vector>
#include <string>

namespace tensorflow
{
    class Session;
    class Tensor;
}

namespace tf
{

class MyGraph
{
public:
   int n_inputs = 1;
   int n_outputs = 1;

   static std::unique_ptr<MyGraph> Mycreate(const char* graph_file_name, const std::vector<std::string> & outputs = {}, int ninputs = 1, int noutputs = 1)
    {
        bool success;
        std::unique_ptr<MyGraph> ptr(new MyGraph(graph_file_name, outputs, success, ninputs, noutputs));
        if (success) { return ptr; }
        else { return nullptr; }
    }

    ~MyGraph();

    std::vector<float> Myrun(const std::vector< std::vector<float> > & x);

    // process vector of 3D inputs, return vector of 1D outputs; use all inputs
    // if samples = -1, or only the specified number of first samples
    std::vector< std::vector < std::vector< float > > > Myrun(
	const std::vector< std::vector< std::vector< std::vector<float> > > > & x,
	long long int samples = -1);
    std::vector< std::vector < std::vector< float > > > Myrun(const std::vector< tensorflow::Tensor > & x);

private:
    /// Not-throwing constructor.
    MyGraph(const char* graph_file_name, const std::vector<std::string> & outputs, bool & success, int ninputs, int noutputs);

    tensorflow::Session* fSession;
    //std::vector< std::string > fInputNames;
    std::vector< std::string > fInputNames;
    std::vector< std::string > fOutputNames;
};

} // namespace tf

#endif
