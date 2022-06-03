cmake_minimum_required(VERSION 3.6)
project(gallery2larcv CXX)

#---Module directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
list(APPEND CMAKE_MODULE_PATH "$ENV{IFDHC_LIB}/ifdhc/cmake")

# dependency hell
find_package(PkgConfig)
pkg_check_modules(SQLite3 REQUIRED sqlite3)
set( SQLite3_LIBRARY $ENV{SQLITE_LIB}/libsqlite3_ups.so )
find_package( Boost REQUIRED )
find_package( ROOT REQUIRED )
set( Range_DIR $CMAKE_CURRENT_SOURCE_DIR )
find_package( Range )
find_package( icarusalg REQUIRED QUIET )
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

add_executable(gallery2larcv gallery2larcv.cxx)

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
include_directories( ${icarusalg_INCLUDE_DIRS} )
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
  icarusalg::Geometry 
  art::Persistency_Provenance
  art::Persistency_Common
  art::Framework_Principal
  nusimdata::SimulationBase
  lardataobj::RecoBase
  lardataobj::RawData
  lardataobj::AnalysisBase
  lardataobj::OpticalDetectorData
  lardataobj::MCBase
  lardataalg::DetectorInfo
)
set( LARLITE_LIBS
     larlite_LArUtil
     larlite_Base
     larlite_DataFormat
     litemaker_scanner )

set( LINK_LIBS ${ART_LIBS}  ${LARLITE_LIBS} LArCVCoreDataFormat ${ROOT_LIBRARIES} )

#target_link_libraries(demo ${LINK_LIBS} )
target_link_libraries(gallery2larcv ${LINK_LIBS})

