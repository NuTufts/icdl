
set( LARLITE_LIB -L$ENV{LARLITE_LIBDIR} -llarlite_Base -llarlite_DataFormat -llarlite_LArUtil)

include_directories ( $ENV{LARLITE_BASEDIR} )

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} ${LARLITE_LIB}) # All execs need this.

art_make( EXCLUDE LLMetaMaker_service.cc UtilScanner_module.cc
	  BASENAME_ONLY
	  LIBRARY_NAME litemaker_LiteMaker
	  LIB_LIBRARIES
	  	   larcorealg_Geometry
		   larcore_Geometry_Geometry_service
	   	   larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
		   larsim_MCCheater_ParticleInventoryService_service
		   lardataobj_RecoBase
		   nusimdata_SimulationBase
		   lardataobj_MCBase
		   larcoreobj_SummaryData
		   lardataobj_RawData
		   lardata_RecoObjects
		   nusimdata_SimulationBase
		   ${LARRECO_LIB}
		   ${LARDATA_LIB}
		   ${ART_FRAMEWORK_CORE}
		   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_BASIC}
		   ${ART_FRAMEWORK_SERVICES_REGISTRY}
		   ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		   ${ART_ROOT_IO_TFILESERVICE_SERVICE}
		   art_Persistency_Common canvas
		   art_Persistency_Provenance canvas
		   art_Utilities canvas
		   ${MF_MESSAGELOGGER}
		   
		   ${FHICLCPP}
		   cetlib cetlib_except
		   ${ROOT_GEOM}
		   ${ROOT_XMLIO}
		   ${ROOT_GDML}
		   ${ROOT_BASIC_LIB_LIST}

          MODULE_LIBRARIES 	  		   
		   LLMetaMaker_service
                   lardataalg_DetectorInfo
		   larcorealg_Geometry
		   larcore_Geometry_Geometry_service
		   larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
		   nusimdata_SimulationBase
		   lardata_Utilities
		   lardataobj_MCBase
		   larcoreobj_SummaryData
		   lardataobj_RawData
		   lardataobj_RecoBase
		   lardata_RecoObjects
                   larevt_Filters
                   larevt_CalibrationDBI_IOVData
                   larevt_CalibrationDBI_Providers
		   nusimdata_SimulationBase
		   ${LARRECO_LIB}
		   ${LARDATA_LIB}
		   ${ART_FRAMEWORK_CORE}
		   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_BASIC}
		   ${ART_FRAMEWORK_SERVICES_REGISTRY}
		   ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		   ${ART_ROOT_IO_TFILESERVICE_SERVICE}
		   art_Persistency_Common canvas
		   art_Persistency_Provenance canvas
		   art_Utilities canvas
		   ${MF_MESSAGELOGGER}
		   
		   ${FHICLCPP}
		   cetlib cetlib_except
		   ${ROOT_GEOM}
		   ${ROOT_XMLIO}
		   ${ROOT_GDML}
		   ${ROOT_BASIC_LIB_LIST}
		)

simple_plugin(LLMetaMaker "service"
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}	      
              art_Framework_Principal
	      art_Utilities canvas
	      ${MF_MESSAGELOGGER}
	      
	      ${FHICLCPP}
	      cetlib cetlib_except
	      ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

# install_headers()
#install_fhicl()
#install_source()
