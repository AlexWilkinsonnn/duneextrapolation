add_subdirectory(run_fcls)

art_make(
        MODULE_LIBRARIES

                lardataobj_RecoBase
                lardataobj_RawData
                lardataobj_Simulation
                nusimdata_SimulationBase
                larcorealg_Geometry
                larcore_Geometry_Geometry_service
                lardataalg_DetectorInfo
                larreco_Calorimetry

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

                ${ROOT_GEOM}
                ${ROOT_XMLIO}
                ${ROOT_GDML}
                ${ROOT_BASIC_LIB_LIST}
                ${ROOT_SPECTRUM}
                ${ROOT_ROOFIT}
                ${ROOT_ROOFITCORE}

                dunereco_CVN_func

                hep_hpc_hdf5
                ${HDF5_LIBRARIES}
        )

install_headers()
install_fhicl()
install_source()
