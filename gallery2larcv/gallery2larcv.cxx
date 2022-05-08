/**
 * @file    gallery2larcv.cpp
 * @brief   Simple conversion of ICARUS larsoft file into DL formats (larcv, larlite)
 * @date    May 8, 2022
 * 
 * 
 */

// - Geometry
#include "larcorealg/Geometry/StandaloneGeometrySetup.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "icarusalg/Geometry/ICARUSChannelMapAlg.h"
#include "icarusalg/Geometry/LoadStandardICARUSgeometry.h"
// - configuration
#include "larcorealg/Geometry/StandaloneBasicSetup.h"

// // gallery/canvas
#include "fhiclcpp/ParameterSet.h"

// ROOT
#include "TFile.h"
#include "TTree.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <memory> // std::make_unique()
#include <iostream> // std::cerr
#include <algorithm> // std::copy()
#include <iterator> // std::back_inserter()

#if !defined(__CLING__)
int main(int argc, char** argv) {
  
  char **pParam = argv + 1, **pend = argv + argc;
  if (pParam == pend) {
    std::cerr << "Usage: " << argv[0] << "  configFile"
      << std::endl;
    return 1;
  }
  std::string const configFile = *(pParam++);
  std::vector<std::string> fileNames;
  std::copy(pParam, pend, std::back_inserter(fileNames));

  // read FHiCL configuration from a configuration file:
  fhicl::ParameterSet config = lar::standalone::ParseConfiguration(configFile);
  
  // // set up message facility (always picked from "services.message")
  lar::standalone::SetupMessageFacility(config, "servicesdemo");
    
  //geometry setup
  //auto geom = lar::standalone::SetupGeometry<icarus::ICARUSChannelMapAlg>(config.get<fhicl::ParameterSet>("services.Geometry"));
  //geom->Print(std::cout);
    
  return 0;
} // main()

#endif // !__CLING__

