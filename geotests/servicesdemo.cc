/**
 * @file    galleryAnalysis.cpp
 * @brief   Template analysis program based on gallery.
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    October 21, 2017
 * 
 * The script is an adaptation of the official gallery demo script demo.cc at
 * https://github.com/marcpaterno/gallery-demo .
 * 
 * To jump into the action, look for `SERVICE PROVIDER SETUP` and
 * `SINGLE EVENT PROCESSING` tags in the source code.
 * 
 * The approach for loading services is the lowest level LArSoft provides.
 * An higher level one is to use `testing::TesterEnvironment` as in some service
 * provider unit tests (e.g., `geo::GeometryCore` and `detinfo::LArProperties`).
 * 
 */

// LArSoft
// // - data products
// #include "lardataobj/RecoBase/Track.h"
// // - DetectorProperties
// #include "lardataalg/DetectorInfo/DetectorPropertiesStandardTestHelpers.h"
// #include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
// // - DetectorClocks
// #include "lardataalg/DetectorInfo/DetectorClocksStandardTestHelpers.h"
// #include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
// // - LArProperties
// #include "lardataalg/DetectorInfo/LArPropertiesStandardTestHelpers.h"
// #include "lardataalg/DetectorInfo/LArPropertiesStandard.h"
// - Geometry
#include "larcorealg/Geometry/StandaloneGeometrySetup.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "icarusalg/Geometry/ICARUSChannelMapAlg.h"
#include "icarusalg/Geometry/LoadStandardICARUSgeometry.h"
// - configuration
#include "larcorealg/Geometry/StandaloneBasicSetup.h"

// // gallery/canvas
// #include "gallery/Event.h"
// #include "canvas/Utilities/InputTag.h"
// #include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// // ROOT
// #include "TFile.h"

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

  /*
   * the "test" environment configuration
   */
  // read FHiCL configuration from a configuration file:
  fhicl::ParameterSet config = lar::standalone::ParseConfiguration(configFile);
  
  // // set up message facility (always picked from "services.message")
  lar::standalone::SetupMessageFacility(config, "servicesdemo");
  
  // // configuration from the "analysis" table of the FHiCL configuration file:
  // auto const& analysisConfig = config.get<fhicl::ParameterSet>("analysis");
  
  // // ***************************************************************************
  // // ***  SERVICE PROVIDER SETUP BEGIN  ****************************************
  // // ***************************************************************************
  // //
  // // Uncomment the things you need
  // // (and make sure the corresponding headers are also uncommented)
  // //
  
  //geometry setup (it's special)
  auto geom = lar::standalone::SetupGeometry<icarus::ICARUSChannelMapAlg>(config.get<fhicl::ParameterSet>("services.Geometry"));
  //std::unique_ptr<geo::GeometryCore> geom = icarus::geo::LoadStandardICARUSgeometry("standard_g4_icarus.fcl");
  //geom->Print(std::cout);
  
  // LArProperties setup
  //auto larp = testing::setupProvider<detinfo::LArPropertiesStandard>(config.get<fhicl::ParameterSet>("services.LArPropertiesService"));
  
  // DetectorClocks setup
  //auto detclk = testing::setupProvider<detinfo::DetectorClocksStandard>(config.get<fhicl::ParameterSet>("services.DetectorClocksService"));
  
  // DetectorProperties setup
  // auto detp = testing::setupProvider<detinfo::DetectorPropertiesStandard>(config.get<fhicl::ParameterSet>("services.DetectorPropertiesService"),
  // 									  detinfo::DetectorPropertiesStandard::providers_type{geom.get(),
  // 									    static_cast<detinfo::LArProperties const*>(larp.get()),
  // 									    // TODO type cast is required until issue #18001 is solved
  // 									    static_cast<detinfo::DetectorClocks const*>(detclk.get())});
  
  // ***************************************************************************
  // ***  SERVICE PROVIDER SETUP END    ****************************************
  // ***************************************************************************
  
  
  return 0;
} // main()

#endif // !__CLING__

