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
// - data objects
#include "lardataobj/RecoBase/Wire.h"

// // gallery/canvas
#include "fhiclcpp/ParameterSet.h"
#include "gallery/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

// LARCV
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"

// ROOT
#include "TFile.h"
#include "TTree.h"

// C/C++ standard libraries
#include <cmath>
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
  auto geom = lar::standalone::SetupGeometry<icarus::ICARUSChannelMapAlg>(config.get<fhicl::ParameterSet>("services.Geometry"));

  // get parameters needed for image making
  unsigned int max_wires_per_plane   = geom->MaxWires();
  unsigned int max_ticks_per_channel = 4096;
  std::cout << "max wires per plane: " << max_wires_per_plane << std::endl;  

  std::vector<std::string> allInputFiles
    //= { "../../testdata/mcc2021b_bnbnu/prodcorsika_bnb_genie_protononly_overburden_icarus_gen_filter_g4_detsim_24027129_667_reco1.root"};
    = {"/icarus/data/users/tmw/first_mc_files/prodcorsika_bnb_genie_protononly_overburden_icarus_gen_filter_g4_detsim_24027129_667_reco1.root"};

  // OUTPUT SETUP
  larcv::IOManager iolcv( larcv::IOManager::kWRITE );
  iolcv.set_out_file( "out_larcv.root" );
  iolcv.initialize();

  int numEvents(0);

  art::InputTag const wire_tag("roifinder","PHYSCRATEDATATPCEW","MCstage0");
  
  std::cout << "start event loop" << std::endl;

  /*
   * the event loop
   */
  int ientry = -1;
  for (gallery::Event event(allInputFiles); !event.atEnd(); event.next())
    {
      ientry++;
      
      // *************************************************************************
      // ***  SINGLE EVENT PROCESSING BEGIN  *************************************
      // *************************************************************************
      
      mf::LogVerbatim("gallery2larcv") << "This is event " << event.fileEntry() << "-" << event.eventEntry();

      // we get the wire data
      
      auto pWireVec = event.getValidHandle< std::vector<recob::Wire> >(wire_tag);
      int nwires = pWireVec->size();
      std::cout << " num wires in event: " << nwires << std::endl;

      auto ev_wireout = (larcv::EventImage2D*)iolcv.get_data(larcv::kProductImage2D,"wire");
      for (int iplane=0; iplane<24; iplane++) {
	// 2 cryo * 4 tpc * 3 planes = 24 planes total
	larcv::ImageMeta meta( max_wires_per_plane, max_ticks_per_channel, 
			       (int)max_ticks_per_channel, (int)max_wires_per_plane,
			       0.0, 0.0, iplane );
	larcv::Image2D img(meta);
	img.paint(0.0);
	ev_wireout->Emplace( std::move(img) );
      }
      auto& adc_v = ev_wireout->as_mut_vector();

      int nmaxsignal = 0;
      for ( auto const& wire : *pWireVec )  {
	std::vector<float>  signal = wire.Signal();
	raw::ChannelID_t   channel = wire.Channel();
	std::vector<geo::WireID> wire_v = geom->ChannelToWire(channel);
	if ( signal.size()>nmaxsignal )
	  nmaxsignal = signal.size();
	// get plane
	if ( wire_v.size()>0 ) {
	  auto const& wireid = wire_v.at(0);
	  int wid     = wireid.Wire;
	  int planeid = wireid.parentID().Plane;
	  int tpcid  = wireid.parentID().parentID().TPC;
	  int cryoid = wireid.parentID().parentID().parentID().Cryostat;
	  int lcvplaneid = cryoid*(12) + tpcid*3 + planeid;
	  std::cout << "  wire (p,t,c)=(" << planeid << "," << tpcid << "," << cryoid << ")" << std::endl;
	  auto& img = adc_v.at(lcvplaneid);
	  img.copy( 0, wid, signal, std::min((int)signal.size(),(int)4096) );
	}
	//std::cout << "CH " << channel << ": no wires=" << wire_v.size() << std::endl;
      }
      std::cout << "max signal size: " << nmaxsignal << std::endl;
      
      iolcv.set_id( 0, 0, ientry );
      iolcv.save_entry();
      // define a meta
      // we transfer data to image
      // we write output
      if (true)
	break;
    }

    
  //geom->Print(std::cout);
  iolcv.finalize();
    
  return 0;
} // main()

#endif // !__CLING__

