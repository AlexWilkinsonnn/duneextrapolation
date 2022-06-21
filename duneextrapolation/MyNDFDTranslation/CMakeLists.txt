add_subdirectory(run_fcls)

set(CMAKE_PREFIX_PATH "/work1/dune/users/awilkins/extrapolation/pytorch-install")
find_package(Torch REQUIRED)

# using SYSTEM to prevent multiple "error: extra ‘;’ [-Werror=pedantic]" errors
include_directories(SYSTEM ${TORCH_INCLUDE_DIRS})

art_make(
        MODULE_LIBRARIES

                lardataobj_RecoBase
                lardataobj_RawData
                lardataobj_Simulation
                nusimdata_SimulationBase
                larcorealg_Geometry
                larcore_Geometry_Geometry_service
                lardataalg_DetectorInfo
                #dune-raw-data_Overlays
                #dune_ArtSupport

                art_Persistency_Common canvas
                art_Persistency_Provenance canvas
                art_Utilities canvas
                cetlib cetlib_except
                ${LARRECO_LIB}
                ${LARDATA_LIB}
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_FRAMEWORK_SERVICES_BASIC}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                ${MF_MESSAGELOGGER}
                ${MF_UTILITIES}
                ${FHICLCPP}
                ${CLHEP}

                # Complains about absolute paths, not sure If I need this anyway
                # ${TORCH_INSTALL_PREFIX}/lib/libtorch.so
                # ${TORCH_INSTALL_PREFIX}/lib/libtorch_cpu.so
                # ${TORCH_INSTALL_PREFIX}/lib/libtorch_cuda.so
                # ${TORCH_INSTALL_PREFIX}/lib/libc10.so
                # ${TORCH_INSTALL_PREFIX}/lib/libc10_cuda.so
                torch

                ${ROOT_GEOM}
                ${ROOT_XMLIO}
                ${ROOT_GDML}
                ${ROOT_BASIC_LIB_LIST}
                ${ROOT_SPECTRUM}
                ${ROOT_ROOFIT}
                ${ROOT_ROOFITCORE}
        )

install_headers()
install_fhicl()
install_source()
