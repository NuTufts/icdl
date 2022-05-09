set(MODULE_NAME scanner)

set( LIBNAME litemaker_${MODULE_NAME} )

# Collect the headers
set( HEADERS ScannerAlgo.h 
  ScannerAlgo.template.h
)

# Add sources
add_library( ${LIBNAME} SHARED
  ScannerAlgo.cxx
  )

# create a property to hold our headers
set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")
get_target_property(ALL_HEADERS ${LIBNAME} PUBLIC_HEADER)

# includes
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR} )
include_directories( ${fhiclcpp_INCLUDE_DIRS} )
include_directories( ${cetlib_INCLUDE_DIRS} )
include_directories( ${cetlib_except_INCLUDE_DIRS} )
include_directories( ${messagefacility_INCLUDE_DIRS} )
include_directories( ${art_INCLUDE_DIRS} )
include_directories( ${canvas_INCLUDE_DIRS} )
include_directories( ${canvas_root_io_INCLUDE_DIRS} )
include_directories( ${gallery_INCLUDE_DIRS} )
include_directories( ${lardataobj_INCLUDE_DIRS} )
include_directories( ${larcoreobj_INCLUDE_DIRS} )
include_directories( ${larcorealg_INCLUDE_DIRS} )
include_directories( ${nusimdata_INCLUDE_DIRS} )
include_directories( ${LARLITE_INCLUDE_DIR} )

target_include_directories(${LIBNAME}
  PUBLIC
  $<INSTALL_INTERFACE:include/litemaker/${MODULE_NAME}>
  ${fhiclcpp_INCLUDE_DIRS}
  ${cetlib_INCLUDE_DIRS}
  ${cetlib_except_INCLUDE_DIRS}
  ${messagefacility_INCLUDE_DIRS}
  ${art_INCLUDE_DIRS}
  ${canvas_INCLUDE_DIRS}
  ${canvas_root_io_INCLUDE_DIRS}
  ${gallery_INCLUDE_DIRS}
  ${lardataobj_INCLUDE_DIRS}
  ${larcoreobj_INCLUDE_DIRS}
  ${larcorealg_INCLUDE_DIRS}
  ${nusimdata_INCLUDE_DIRS}
  ${LARLITE_INCLUDE_DIR}
  PRIVATE
  ${PROJECT_SOURCE_DIR}
  )


# Generate the dictionary
ROOT_GENERATE_DICTIONARY(G__${MODULE_NAME} ${HEADERS} LINKDEF LinkDef.h)
 
# Generate the shared library from the sources
add_library(${MODULE_NAME} SHARED ${SOURCES} G__${MODULE_NAME}.cxx)
target_link_libraries(${MODULE_NAME} PUBLIC 
  larlite_DataFormat 
  larlite_Base 
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  messagefacility::MF_MessageLogger 
  art::Persistency_Common
  canvas::canvas
  ${ROOT_LIBRARIES})

# install libraries and headers
install(TARGETS ${LIBNAME}
  EXPORT litemaker
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/litemaker/${MODULE_NAME})

# install files needed for dictionary use in CINT/CLING/PYTHON
install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lib${MODULE_NAME}_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/lib${MODULE_NAME}.rootmap
  DESTINATION lib  )
