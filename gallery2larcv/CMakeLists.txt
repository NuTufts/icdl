cmake_minimum_required(VERSION 3.6)
find_package(cetmodules 3.20.00 REQUIRED)

project(gallery2larcv CXX)

#---Target executable we are building
add_executable(gallery2larcv gallery2larcv.cxx)

include(CetCMakeEnv)
cet_cmake_env()

#---Module directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
list(APPEND CMAKE_MODULE_PATH "$ENV{IFDHC_LIB}/ifdhc/cmake")

#--- Are we ICARUS, SBND, UBOONE? ---
if (DEFINED ENV{SBNDCODE_DIR})
  set( EXPERIMENT "SBND" )
  message( STATUS "Experiment: SBND" )
  target_compile_definitions( gallery2larcv PUBLIC SBND )
elseif (DEFINED ENV{ICARUSCODE_DIR})
  set( EXPERIMENT "ICARUS" )
  message( STATUS "Experiment: ICARUS" )
  target_compile_definitions( gallery2larcv PUBLIC ICARUS )
elseif (DEFINED  ENV{UBOONECODE_DIR})
  set( EXPERIMENT "UBOONE" )
  message( STATUS "Experiment: UBOONE" )
  target_compile_definitions( gallery2larcv PUBLIC UBOONE )
else()
  message( FATAL_ERROR "Did not detect any of the experiment codebases have been setup. use UPS 'setup uboonecode/sbndcode/icaruscode'")
endif()
     

# welcome to dependency hell
find_package(PkgConfig)
set( SQLite3_LIBRARY $ENV{SQLITE_LIB}/libsqlite3_ups.so )
pkg_check_modules(SQLite3 REQUIRED sqlite3)
find_package( Boost REQUIRED )
find_package( ROOT REQUIRED )
set( Range_DIR $CMAKE_CURRENT_SOURCE_DIR )
find_package( Range-v3 )

if ( EXPERIMENT STREQUAL "ICARUS" )
   find_package( icarusalg REQUIRED )
elseif (EXPERIMENT STREQUAL "SBND" )
   find_package( sbndcode REQUIRED ) 
endif()

find_package( canvas REQUIRED QUIET )
find_package( gallery REQUIRED QUIET )
find_package( art REQUIRED QUIET )
find_package( cetlib REQUIRED QUIET )
find_package( cetlib_except REQUIRED QUIET )
find_package( lardataalg REQUIRED QUIET )
find_package( lardata REQUIRED QUIET )
find_package( sbnobj REQUIRED QUIET )
set( larlite_DIR $ENV{LARLITE_LIBDIR}/cmake/larlite )
find_package( larlite REQUIRED )
set( larcv_DIR $ENV{LARCV_LIBDIR}/cmake/larcv )
find_package( larcv REQUIRED )
set( litemaker_DIR $ENV{ICDL_BASEDIR}/litemaker/build/installed/lib/cmake/litemaker )
find_package( litemaker REQUIRED )
find_package( ifdhc REQUIRED QUIET )
find_package( nusimdata REQUIRED QUIET )

# To set the C++ standard to be used, only one of the following two lines should
# be un-commented.
#set_property(TARGET demo PROPERTY CXX_STANDARD 14)
set_property(TARGET gallery2larcv PROPERTY CXX_STANDARD 17)

include_directories( ${fhiclcpp_INCLUDE_DIRS} )
include_directories( ${canvas_INCLUDE_DIRS} )
include_directories( ${canvas_root_io_INCLUDE_DIRS} )
include_directories( ${gallery_INCLUDE_DIRS} )
include_directories( ${messagefacility_INCLUDE_DIRS} )
include_directories( ${lardataobj_INCLUDE_DIRS} )
include_directories( ${larcorealg_INCLUDE_DIRS} )
include_directories( ${larcorealg_INCLUDE_DIRS} )
if ( EXPERIMENT} STREQUAL "ICARUS" )
  include_directories( ${icarusalg_INCLUDE_DIRS} )
elseif( EXPERIMENT} STREQUAL "ICARUS" )
  include_directories( ${sbndcode_INCLUDE_DIRS} )
endif()

include_directories( ${LARLITE_INCLUDE_DIR} )
include_directories( ${LARCV_INCLUDE_DIR} )
include_directories( ${LITEMAKER_INCLUDE_DIR} )

include_directories( ${cetlib_INCLUDE_DIRS} )
include_directories( ${cetlib_except_INCLUDE_DIRS} )
include_directories( ${art_INCLUDE_DIRS} )
include_directories( ${nusimdata_INCLUDE_DIRS} )

set( ART_LIBS
  fhiclcpp::fhiclcpp 
  gallery::gallery 
  canvas::canvas 
  larcorealg::Geometry 
  art::Persistency_Provenance
  art::Persistency_Common
  art::Framework_Principal
  nusimdata::SimulationBase
  lardataobj::RecoBase
  lardataobj::RawData
  lardataobj::AnalysisBase
  lardataobj::OpticalDetectorData_dict
  lardataobj::MCBase
  lardataalg::DetectorInfo
)
if ( EXPERIMENT STREQUAL "ICARUS" )
   list( APPEND ART_LIBS icarusalg::Geometry  )
elseif ( EXPERIMENT STREQUAL "SBND" )
   list( APPEND ART_LIBS sbndcode::Geometry  )
endif()

set( LARLITE_LIBS
     larlite_LArUtil
     larlite_Base
     larlite_DataFormat
     litemaker_scanner )

set( LINK_LIBS ${ART_LIBS}  ${LARLITE_LIBS} LArCVCoreDataFormat ${ROOT_LIBRARIES} )

#target_link_libraries(demo ${LINK_LIBS} )
target_link_libraries(gallery2larcv ${LINK_LIBS})

