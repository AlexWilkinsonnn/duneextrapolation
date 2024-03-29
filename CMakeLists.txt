
# ======================================================================
#  duneextrapolation main build file
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

project(duneextrapolation VERSION 09.42.04.00 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

find_ups_boost( )

cet_report_compiler_flags()

find_ups_product( eigen )
include_directories( $ENV{EIGEN_INC} )
find_ups_product( art )
find_ups_product( art_root_io )
find_ups_product( canvas_root_io )
find_ups_boost( )
find_ups_root()
find_ups_product( cetbuildtools )
find_ups_product( dunecore )
find_ups_product( dunereco )
find_ups_product( duneopdet )
find_ups_product( lardataalg )
find_ups_product( lardataobj )
find_ups_product( larcore )
find_ups_product( larcorealg )
find_ups_product( larcoreobj )
find_ups_product( nurandom )
find_ups_product( larreco )
find_ups_product( larsim )
find_ups_product( larevt )
find_ups_product( lardata )
find_ups_product( larpandora )
find_ups_product( systematicstools )
find_ups_product( larrecodnn )
find_ups_product( larpandora )
find_ups_product( clhep )
find_ups_geant4( )
if(DEFINED ENV{CAFFE_LIB} )
  find_ups_product(caffe)
endif()
if(DEFINED ENV{LAPACK_LIB} )
  find_ups_product( lapack )
endif()
if(( DEFINED ENV{TENSORFLOW_DIR} ) AND ( DEFINED ENV{TRITON_DIR} ))
  # find tensorflow library
  find_ups_product( protobuf )
  find_ups_product( tensorflow )
  find_ups_product( triton )
  find_ups_product( larrecodnn )
  # cet_find_library(TENSORFLOW_CC NAMES tensorflow_cc PATHS ENV TENSORFLOW_LIB NO_DEFAULT_PATH )
  # cet_find_library(TENSORFLOW_FRAMEWORK NAMES tensorflow_framework PATHS ENV TENSORFLOW_LIB NO_DEFAULT_PATH )
  # set(TENSORFLOW TENSORFLOW_CC TENSORFLOW_FRAMEWORK)
  # cet_find_library(PROTOBUF NAMES protobuf PATHS ENV PROTOBUF_LIB NO_DEFAULT_PATH )
endif()
find_ups_product( hdf5 )
if(DEFINED ENV{HEP_HPC_DIR} )
  find_ups_product( hep_hpc )
endif()
find_package( HighFive )
include_directories( $ENV{HIGHFIVE_INC} )

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(duneextrapolation)
add_subdirectory(ups)
add_subdirectory(fcl)

include(UseCPack)
