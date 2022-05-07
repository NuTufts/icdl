/**
 * @file    galleryAnalysis.cpp
 * @brief   Template analysis program based on gallery.
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    October 21, 2017
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

#include "larlite/Base/GeoTypes.h"
#include "larlite/Geometry/CryoGeo.h"

larlite::geo::View_t convert_view( geo::View_t v)
{
  switch(v) {
  case geo::kU:
    return larlite::geo::kU;
  case geo::kV:
    return larlite::geo::kV;
  case geo::kW:
    return larlite::geo::kW;
  case geo::kY:
    return larlite::geo::kY;
  case geo::kX:
    return larlite::geo::kX;
  case geo::k3D:
    return larlite::geo::k3D;
  default:
    return larlite::geo::kUnknown;
  }//end of switch
  return larlite::geo::kUnknown;
}

larlite::geo::SigType_t convert_signal( geo::SigType_t t ) 
{
  switch (t) {
  case geo::kInduction:
    return larlite::geo::kInduction;
  case geo::kCollection:
    return larlite::geo::kCollection;
  default:
    return larlite::geo::kMysteryType;
    break;
  }
  return larlite::geo::kMysteryType;
}

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

  TFile outfile("icarus_larlite_geodata.root","recreate");
  TTree geodata("geodata","Geometry data for ICARUS detector");
  std::vector< larlite::geo::CryoGeo > cryo_v;
  geodata.Branch( "cryo_v", &cryo_v );
    
  //geometry setup
  auto geom = lar::standalone::SetupGeometry<icarus::ICARUSChannelMapAlg>(config.get<fhicl::ParameterSet>("services.Geometry"));
  //geom->Print(std::cout);
  
  std::cout << "Detector " << geom->DetectorName() << " has "
	    << geom->Ncryostats() << " cryostats and "
	    << geom->NAuxDets() << " auxiliary detectors:";

  // auto const& detEnclosureBox = DetectorEnclosureBox();
  // out << "\n" << indent << "Detector enclosure: "
  //   << detEnclosureBox.Min() << " -- " << detEnclosureBox.Max()
  //   << " cm => ( " << detEnclosureBox.SizeX() << " x "
  //   << detEnclosureBox.SizeY() << " x "
  //   << detEnclosureBox.SizeZ() << " ) cm^3"
  //   ;

  std::string indent = "  ";

  for (geo::CryostatGeo const& cryostat: geom->IterateCryostats()) {
    std::cout << "\n" << indent;
    cryostat.PrintCryostatInfo(std::cout, indent + "  ", cryostat.MaxVerbosity);

    larlite::geo::CryoGeo llcryo( (int)cryostat.ID().Cryostat );
    auto const& cryocenter = cryostat.GetCenter();
    llcryo.fCenter = TVector3( cryocenter.x(), cryocenter.y(), cryocenter.z() );
    auto const& cryomin = cryostat.BoundingBox().Min();
    auto const& cryomax = cryostat.BoundingBox().Max();
    llcryo.fBounds.push_back( TVector3(cryomin.x(),cryomin.y(),cryomin.z()) );
    llcryo.fBounds.push_back( TVector3(cryomax.x(),cryomax.y(),cryomax.z()) );
    
    const unsigned int nTPCs = cryostat.NTPC();
    for(unsigned int t = 0;  t < nTPCs; ++t) {
      const geo::TPCGeo& tpc = cryostat.TPC(t);

      std::cout << "\n" << indent << "  ";
      tpc.PrintTPCInfo(std::cout, indent + "    ", tpc.MaxVerbosity);

      larlite::geo::TPCGeo lltpc( t, (int)cryostat.ID().Cryostat );
      lltpc.fCenter = tpc.GetCenter();
      lltpc.fHalfLengths(0) = tpc.ActiveHalfWidth();
      lltpc.fHalfLengths(1) = tpc.ActiveHalfHeight();
      lltpc.fHalfLengths(2) = tpc.ActiveHalfLength();
      lltpc.fDriftDirection = tpc.DriftDir();
      auto const& tpcmin = tpc.BoundingBox().Min();
      auto const& tpcmax = tpc.BoundingBox().Max();
      lltpc.fBounds.push_back( TVector3(tpcmin.x(),tpcmin.y(),tpcmin.z()) );
      lltpc.fBounds.push_back( TVector3(tpcmax.x(),tpcmax.y(),tpcmax.z()) );
      
      const unsigned int nPlanes = tpc.Nplanes();
      for(unsigned int p = 0; p < nPlanes; ++p) {
        const geo::PlaneGeo& plane = tpc.Plane(p);
        const unsigned int nWires = plane.Nwires();
	
        std::cout << "\n" << indent << "    ";
        plane.PrintPlaneInfo(std::cout, indent + "      ", plane.MaxVerbosity);
        geo::SigType_t const sigType = geom->SignalType(plane.ID());
        std::cout << "\n" << indent << "      "
		  << "signal type: " << SignalTypeName(sigType)
		  << " (" << static_cast<int>(sigType) << ")";
	
	larlite::geo::PlaneGeo llplane( p, t, (int)cryostat.ID().Cryostat, 
					convert_view(plane.View()), 
					convert_signal(geom->SignalType(plane.ID())) );
	llplane.fCenter        = plane.GetCenter();
	llplane.fNormToCathode = plane.GetNormalDirection();
	llplane.fWirePitchDir  = plane.GetIncreasingWireDirection();
	llplane.fWirePitchLen  = plane.WirePitch();
	auto const& planemin = plane.BoundingBox().Min();
	auto const& planemax = plane.BoundingBox().Max();
	llplane.fBoundingBox.push_back( TVector3(planemin.x(),planemin.y(),planemin.z()) );
	llplane.fBoundingBox.push_back( TVector3(planemax.x(),planemax.y(),planemax.z()) );

        for(unsigned int w = 0;  w < nWires; ++w) {
          const geo::WireGeo& wire = plane.Wire(w);
          geo::WireID wireID(plane.ID(), w);
	  
	  raw::ChannelID_t ch_id = geom->PlaneWireToChannel( wireID );
	  
          // the wire should be aligned on z axis, half on each side of 0,
          // in its local frame
          std::cout << "\n" << indent << "      " << wireID << " CH=" << ch_id << " ";
          wire.PrintWireInfo(std::cout, indent + "      ", wire.MaxVerbosity);

	  larlite::geo::WireGeo llwire( (int)ch_id, (int)wireID.Wire, p, t, (int)cryostat.ID().Cryostat,
					wire.GetStart(), wire.GetEnd() );

	  std::cout << indent << "Making larlite::geo::WireGeo[ " << (int)ch_id << "," 
		    <<  (int)wireID.Wire << ","
		    <<  p << "," <<  t << "," <<  (int)cryostat.ID().Cryostat << "]" 
		    << std::endl;

	  llplane.fWires_v.emplace_back( std::move(llwire) );

        } // for wire
	lltpc.planes_v.emplace_back( std::move(llplane) );
      } // for plane in a tpc

      llcryo.tpc_v.emplace_back( std::move(lltpc) );
    } // for TPC
    
    // unsigned int nOpDets = cryostat.NOpDet();
    // for (unsigned int iOpDet = 0; iOpDet < nOpDets; ++iOpDet) {
    //   geo::OpDetGeo const& opDet = cryostat.OpDet(iOpDet);
    //   std::cout << "\n" << indent << "  [OpDet #" << iOpDet << "] ";
    //   opDet.PrintOpDetInfo
    //     (std::forward<Stream>(std::cout), indent + "  ", opDet.MaxVerbosity);
    // } // for

    cryo_v.push_back( std::move(llcryo) );
  } // for cryostat
  
  geodata.Fill();

  std::cout << '\n';
  outfile.Write();
  
  return 0;
} // main()

#endif // !__CLING__

