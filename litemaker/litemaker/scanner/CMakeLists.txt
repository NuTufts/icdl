set(MODULE_NAME scanner)

set( LIBNAME larlite_${MODULE_NAME} )

# Collect the headers
set( HEADERS ScannerAlgo.h
  )

# Add sources
add_library( ${LIBNAME} SHARED
  ScannerAlgo.cxx
  )

# Remove LinkDef.h
list(FILTER HEADERS EXCLUDE REGEX ".*LinkDef.h$")
 
# includes
include_directories( ${LARLITE_INCLUDE_DIRS} )

# Generate the dictionary
ROOT_GENERATE_DICTIONARY(G__${MODULE_NAME} ${HEADERS} LINKDEF LinkDef.h)
 
# Generate the shared library from the sources
add_library(${MODULE_NAME} SHARED ${SOURCES} G__${MODULE_NAME}.cxx)
set_target_properties(${MODULE_NAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_link_libraries(${MODULE_NAME} PUBLIC larlite_DataFormat ${ROOT_LIBRARIES})

# install libraries and headers
install(TARGETS ${MODULE_NAME}
  EXPORT litemaker
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/litemaker/scanner)

# install files needed for dictionary use in CINT/CLING/PYTHON
install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lib${MODULE_NAME}_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/lib${MODULE_NAME}.rootmap
  DESTINATION lib  )
