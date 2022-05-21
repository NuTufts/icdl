/**
 * @file    galleryAnalysis.cpp
 * @brief   Template analysis program based on gallery.
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    October 21, 2017
 * 
 * 
 */


// // gallery/canvas
#include "fhiclcpp/ParameterSet.h"

// larsoft
#include "larcorealg/Geometry/StandaloneBasicSetup.h"

#ifdef ICARUS
#include "icarusalg/Geometry/ICARUSChannelMapAlg.h"
#include "icarusalg/Geometry/LoadStandardICARUSgeometry.h"
#endif

#ifdef SBND
#include "sbndcode/Geometry/ChannelMapSBNDAlg.h"
#endif


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

// litemaker
#include "litemaker/scanner/UtilScannerAlgo.h"


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

  std::cout << "LOAD GEOMETRY" << std::endl;

#ifdef ICARUS
  auto geom = icarus::geo::LoadStandardICARUSgeometry( configFile );
#endif

#ifdef ICARUS
  TFile outfile("icarus_larlite_servicesdata.root","recreate");
#endif
#ifdef SBND
  TFile outfile("sbnd_larlite_servicesdata.root","recreate");
#endif

  litemaker::scanner::UtilScannerAlgo scanner;
  // Get Parameters
  auto larp_pset   = config.get<fhicl::ParameterSet>("services.LArPropertiesService");
  auto lardet_pset = config.get<fhicl::ParameterSet>("services.DetectorClocksService");
  
    std::cout << "RUN SCANNER" << std::endl;
  scanner.SaveDetectorProperties( config, *geom );
  scanner.SaveLArProperties( config, *geom );  

  std::cout << '\n';
  outfile.Write();
  
  return 0;
} // main()

#endif // !__CLING__

