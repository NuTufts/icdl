#ifndef __LITEMAKER_SCANNER_UTILSCANNERALGO_H__
#define __LITEMAKER_SCANNER_UTILSCANNERALGO_H__

// art/gallery/canvas
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"

namespace litemaker {
namespace scanner {

  class UtilScannerAlgo {

  public:

    UtilScannerAlgo() {};
    virtual ~UtilScannerAlgo() {};

    void SaveDetectorProperties(const fhicl::ParameterSet& config,
				const ::geo::GeometryCore& geom );
    

    
  };
  
}
}


#endif
