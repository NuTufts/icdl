cmake_minimum_required(VERSION 3.6)

project(gallerydemo CXX)

set( DETECTOR "ICARUS" )
#set( DETECTOR "SBND" )
message("GEO/SERVICES tests")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

find_package(cetbuildtools REQUIRED)
include(CetCMakeEnv)

# dependency hell
find_package(Threads)
find_package(PkgConfig)
pkg_check_modules(SQLite3 REQUIRED sqlite3)
set( SQLite3_LIBRARY $ENV{SQLITE_LIB}/libsqlite3_ups.so )
find_package( Boost REQUIRED )
find_package( art REQUIRED QUIET )
find_package( canvas REQUIRED QUIET )
find_package( ROOT REQUIRED )
set( Range_DIR $CMAKE_CURRENT_SOURCE_DIR )
find_package( Range )
set( larlite_DIR $ENV{LARLITE_LIBDIR}/cmake/larlite )
find_package( lardataobj REQUIRED )
find_package( lardataalg REQUIRED )
find_package( larcore REQUIRED )
find_package( larlite REQUIRED )
set( litemaker_DIR $ENV{ICDL_BASEDIR}/litemaker/build/installed/lib/cmake/litemaker )
find_package( litemaker REQUIRED )

if ( ${DETECTOR} STREQUAL "SBND" )
message( "BUILD FOR SBND" )
add_compile_definitions( SBND )
find_package( sbndcode REQUIRED QUIET )
endif()

if ( ${DETECTOR} STREQUAL "ICARUS" )
message( "BUILD FOR ICARUS" )
add_compile_definitions( ICARUS )
find_package( icarusalg REQUIRED QUIET )
endif()


#add_executable(demo demo.cc analyze.cc)
add_executable(servicesdemo servicesdemo.cc)
add_executable(gen_geodata_rootfile gen_geodata_rootfile.cc)
add_executable(gen_services_rootfile gen_services_rootfile.cc)

# To set the C++ standard to be used, only one of the following two lines should
# be un-commented.
#set_property(TARGET demo PROPERTY CXX_STANDARD 14)
set_property(TARGET servicesdemo PROPERTY CXX_STANDARD 17)
set_property(TARGET gen_geodata_rootfile PROPERTY CXX_STANDARD 17)
set_property(TARGET gen_services_rootfile PROPERTY CXX_STANDARD 17)

include_directories( ${fhiclcpp_INCLUDE_DIRS} )
include_directories( ${messagefacility_INCLUDE_DIRS} )
include_directories( ${canvas_INCLUDE_DIRS} )
include_directories( ${lardataobj_INCLUDE_DIRS} )
include_directories( ${lardataalg_INCLUDE_DIRS} )
include_directories( ${larcorealg_INCLUDE_DIRS} )
include_directories( ${larcore_INCLUDE_DIRS} )
include_directories( ${icarusalg_INCLUDE_DIRS} )
include_directories( ${LARLITE_INCLUDE_DIR} )
include_directories( ${LITEMAKER_INCLUDE_DIR} )

message( "larlite include dirs: ${LARLITE_INCLUDE_DIR}" )

set( LARLITE_LIBS
     larlite_LArUtil
     larlite_Base
     larlite_DataFormat
     litemaker_scanner )

set( SERVICES_ART_LIBS
     fhiclcpp::fhiclcpp 
     larcorealg::Geometry 
     lardataobj::RecoBase
     lardataobj::MCBase
     lardataalg::DetectorInfo
     #lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
     #lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
     #lardata::DetectorInfoServices_LArPropertiesServiceStandard_service
     canvas::canvas
     art::Utilities
     art::Framework_Principal
     art::Persistency_Provenance
     art::Persistency_Common )


if (DETECTOR STREQUAL "SBND")
set( LINK_LIBS fhiclcpp::fhiclcpp larcorealg::Geometry sbndcode::Geometry larlite_LArUtil larlite_Base ${ROOT_LIBRARIES} )
endif()

if (DETECTOR STREQUAL "ICARUS")
set( LINK_LIBS fhiclcpp::fhiclcpp larcorealg::Geometry icarusalg::Geometry larlite_LArUtil larlite_Base ${ROOT_LIBRARIES} )
endif()

#target_link_libraries(demo ${LINK_LIBS} )
target_link_libraries(servicesdemo ${SERVICES_ART_LIBS} ${LINK_LIBS} )
target_link_libraries(gen_geodata_rootfile ${LINK_LIBS})
target_link_libraries(gen_services_rootfile ${LARLITE_LIBS} ${SERVICES_ART_LIBS} ${LINK_LIBS})

