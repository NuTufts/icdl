/**
 * @file    gallery2larcv.cpp
 * @brief   Simple conversion of SBND larsoft file into DL formats (larcv, larlite)
 * @date    May 8, 2022
 * 
 * 
 */

// - Geometry
#include "larcorealg/Geometry/StandaloneGeometrySetup.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "sbndcode/Geometry/ChannelMapSBNDAlg.h"

// - configuration
#include "larcorealg/Geometry/StandaloneBasicSetup.h"
// - data objects
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/SimChannel.h"

// // gallery/canvas
#include "fhiclcpp/ParameterSet.h"
#include "gallery/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

// LARCV
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"

// larlite
#include "larlite/LArUtil/LArUtilConfig.h"
#include "larlite/LArUtil/Geometry.h"
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/mcpart.h"
#include "larlite/DataFormat/mctruth.h"
#include "larlite/DataFormat/mcflux.h"
#include "larlite/DataFormat/gtruth.h"
#include "larlite/DataFormat/opflash.h"
#include "larlite/DataFormat/crttrack.h"
#include "larlite/DataFormat/crthit.h"

// litemaker
#include "litemaker/scanner/ScannerAlgo.h"

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

  bool HAS_MC  = true;
  bool HAS_CRT = false;
  bool HAS_OPFLASH = false;
  bool HAS_RECO1 = false;
  int  downsample_factor =  4;

  // read FHiCL configuration from a configuration file:
  fhicl::ParameterSet config = lar::standalone::ParseConfiguration(configFile);
  
  // // set up message facility (always picked from "services.message")
  lar::standalone::SetupMessageFacility(config, "servicesdemo");

  //geometry setup
  std::cout << "Load SBND Geometry" << std::endl;
  auto geom = lar::standalone::SetupGeometry<geo::ChannelMapSBNDAlg>(config.get<fhicl::ParameterSet>("services.Geometry"));
  larutil::LArUtilConfig::SetDetector( larlite::geo::kSBND );

  auto llgeom = larlite::larutil::Geometry::GetME();
  std::cout << "larlite geometry loaded" << std::endl;

  const int nCryos = llgeom->Ncryostats();
  int nTPCs   = 0;
  int nPlanes = 0;
  for (int icryo=0; icryo<nCryos; icryo++) {
    for ( int itpc=0; itpc<(int)llgeom->NTPCs(icryo); itpc++ ) {
      nTPCs++;
      for (int iplane=0; iplane<(int)llgeom->Nplanes(itpc,icryo); iplane++) {
	nPlanes++;
      }
      std::cout << "Cryo[" << icryo << "] TPC[" << itpc << "] nplanes=" << llgeom->Nplanes(itpc,icryo) << std::endl;
    }
  }
  std::cout << "Total number of cryostats: " << nCryos << std::endl;
  std::cout << "Total number of TPCs: " << nTPCs << std::endl;
  std::cout << "Total number of Planes: " << nPlanes << std::endl;

  // ==== LARSOFT PRODUCTS ===============
  // Make a list of larsoft products we want to translate into larlite

  // Deconvolved wire signals
  std::string wireproducer = "Reco1";
  std::string stage_name = "caldata";
  std::vector< std::string > wire_tagname_v = {""};

  std::vector< art::InputTag > wire_tag_v;
  for (auto& tagname : wire_tagname_v ) {
    art::InputTag wire_tag( stage_name, tagname, wireproducer );
    wire_tag_v.push_back( wire_tag );
  }

  // generator stage
  art::InputTag const mctruth_genie_tag("generator","","GenieGen");
  art::InputTag const mctruth_corsika_tag("corsika","","GenieGen");
  art::InputTag const mcflux_tag("generator","","GenieGen");
  art::InputTag const gtruth_tag("generator","","GenieGen");
  // g4 stage
  art::InputTag const simch_tag("simdrift","","G4");
  art::InputTag const mcreco_track_tag("mcreco","","G4");
  art::InputTag const mcreco_shower_tag("mcreco","","G4");
  // DetSim
  art::InputTag const opdaq_tag("opdaq","","DetSim");
  // reco 1
  art::InputTag const gaushit_tag("gaushit","","Reco1");
  

  // ==== LARLITE SCANNER ====================
  // litemaker scanner: converts larsoft data products into larlite ones

  larlite::ScannerAlgo scanner; 

  // generator
  scanner.Register( "generator", larlite::data::kMCTruth );
  scanner.Register( "generator", larlite::data::kMCFlux );
  scanner.Register( "generator", larlite::data::kGTruth );
  // G4
  scanner.Register( "largeant",  larlite::data::kMCParticle );
  scanner.Register( "simdrift",  larlite::data::kSimChannel );
  scanner.Register( "mcreco",    larlite::data::kMCTrack );
  scanner.Register( "mcreco",    larlite::data::kMCShower );
  // Detsim
  scanner.Register( "opdaq",     larlite::data::kOpDetWaveform );
  // Reco 1
  scanner.Register( "gaushit" ,  larlite::data::kHit );

  larlite::storage_manager llout( larlite::storage_manager::kWRITE );
  llout.set_out_filename( "out_larlite.root" );
  llout.open();

  // get parameters needed for image making
  unsigned int max_wires_per_plane   = geom->MaxWires();
  std::cout << "max wires per plane: " << max_wires_per_plane << std::endl;  

  std::vector<std::string> allInputFiles
    = {"/sbnd/data/users/tmw/test_data/2022A_rockbox/gen_g4_detsim_reco1-69432ad3-6cd6-4591-be56-7fc0adb2c2f1.root"};

  // OUTPUT SETUP
  larcv::IOManager iolcv( larcv::IOManager::kWRITE );
  iolcv.set_out_file( "out_larcv.root" );
  iolcv.initialize();

  int numEvents(0);

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
      
      mf::LogVerbatim("gallery2larcv_sbnd") << "This is event " << event.fileEntry() << "-" << event.eventEntry();

      auto ev_wireout = (larcv::EventImage2D*)iolcv.get_data(larcv::kProductImage2D,"wire");

      std::vector< larcv::Image2D > adc_v;

      // we get the wire data, find longest signal size
      int max_ticks_per_channel = 0;
      for ( auto& wire_tag : wire_tag_v ) {
	auto pWireVec = event.getValidHandle< std::vector<recob::Wire> >(wire_tag);
	int nwires = pWireVec->size();

	for ( auto const& wire : *pWireVec )  {
	  size_t nsignal = wire.Signal().size();
	  if ( (int)nsignal>max_ticks_per_channel )
	    max_ticks_per_channel = (int)nsignal;
	}

      }
      std::cout << "Max ticks per channel: " << max_ticks_per_channel << std::endl;
	  
      int output_nticks = max_ticks_per_channel/downsample_factor;
      int output_nticks_pre = max_ticks_per_channel;
      if ( downsample_factor>1 && max_ticks_per_channel%downsample_factor!=0 ) {
	output_nticks += 1;
	output_nticks_pre = downsample_factor*output_nticks;
      }
      std::cout << "Output ticks per channel, fit with downsampler : " << output_nticks << std::endl;
      std::cout << "Output ticks per channel, pre-downsampler : " << output_nticks_pre << std::endl;

      // we allocate the space for the images.
      for (int icryo=0; icryo<(int)llgeom->Ncryostats(); icryo++ ) {
	for (int itpc=0; itpc<(int)llgeom->NTPCs(icryo); itpc++) {
	  for (int iplane=0; iplane<(int)llgeom->Nplanes(itpc,icryo); iplane++) {
	    int index = llgeom->GetSimplePlaneIndexFromCTP( icryo, itpc, iplane );

	    larcv::ImageMeta meta( max_wires_per_plane, output_nticks_pre,
				   (int)output_nticks_pre, (int)max_wires_per_plane,
				   0.0, 0.0, index );
	    larcv::Image2D img(meta);
	    img.paint(0.0);
	    std::cout << "blank made: " << meta.dump() << std::endl;
	    adc_v.emplace_back( std::move(img) );
	  }
	}
      }
      
      // populate the images
      for ( auto& wire_tag : wire_tag_v ) {

	auto pWireVec = event.getValidHandle< std::vector<recob::Wire> >(wire_tag);
	int nwires = pWireVec->size();

	for ( auto const& wire : *pWireVec )  {
	  std::vector<float>  signal = wire.Signal();
	  raw::ChannelID_t   channel = wire.Channel();
	  std::vector<geo::WireID> wire_v = geom->ChannelToWire(channel);

	  if ( signal.size()==0 )
	    continue;

	  // get plane
	  if ( wire_v.size()>0 ) {
	    auto const& wireid = wire_v.at(0);
	    int wid     = wireid.Wire;
	    int planeid = wireid.parentID().Plane;
	    int tpcid  = wireid.parentID().parentID().TPC;
	    int cryoid = wireid.parentID().parentID().parentID().Cryostat;
	    int lcvplaneid = llgeom->GetSimplePlaneIndexFromCTP( cryoid, tpcid, planeid );
	    
	    //std::cout << "  wire (p,t,c)=(" << planeid << "," << tpcid << "," << cryoid << ")" << std::endl;
	    auto& img = adc_v.at(lcvplaneid);
	    //if ( signal.size()>max_ticks_per_channel ) {
	    //std::cout << "CH " << channel << ": num wires=" << wire_v.size() << " num ticks=" << signal.size() << std::endl;
	    //}
	    img.copy( 0, wid, signal, std::min((int)signal.size(),(int)output_nticks_pre) );
	    //std::cout << "CH " << channel << ": no wires=" << wire_v.size() << std::endl;
	  }
	  //std::cout << "max signal size: " << nmaxsignal << std::endl;
	}//end of loop over wire-vector
      }//end of wire-tag loop
	
      for (auto& img : adc_v ) {

	// compress image first
	if ( downsample_factor>1  )
	  img.compress( output_nticks, max_wires_per_plane, larcv::Image2D::kSum );
	
	int nabove = 0;
	for ( auto const& pixval : img.as_vector() ) {
	  if ( pixval>5.0 )
	    nabove++;
	  
	  if (nabove>10) {
	    break;
	  }
	}
	  
	if ( nabove>10 ) {
	  std::cout << "saving image with content: " << img.meta().dump() << std::endl;
	  ev_wireout->Emplace( std::move(img) );
	}
      }
      
      iolcv.set_id( 0, 0, ientry );
      iolcv.save_entry();
      
      // larlite products
      if ( HAS_MC ) {

	auto pMCTruthVec = event.getValidHandle< std::vector<simb::MCTruth> >(mctruth_genie_tag);
	larlite::event_mctruth* ev_mctruth = 
	  (larlite::event_mctruth*)llout.get_data( larlite::data::kMCTruth, "generator" );
	scanner.ScanData( *pMCTruthVec, ev_mctruth );

	auto pMCFluxVec = event.getValidHandle< std::vector<simb::MCFlux> >(mcflux_tag);
	larlite::event_mcflux* ev_mcflux = 
	  (larlite::event_mcflux*)llout.get_data( larlite::data::kMCFlux, "generator" );
	scanner.ScanData( *pMCFluxVec, ev_mcflux );

	auto pGTruthVec = event.getValidHandle< std::vector<simb::GTruth> >(gtruth_tag);
	larlite::event_gtruth* ev_gtruth = 
	  (larlite::event_gtruth*)llout.get_data( larlite::data::kGTruth, "generator" );
	scanner.ScanData( *pGTruthVec, ev_gtruth );

	auto pMCRecoTrack = event.getValidHandle< std::vector<sim::MCTrack> >( mcreco_track_tag );
	larlite::event_mctrack* ev_mctrack =
	  (larlite::event_mctrack*)llout.get_data( larlite::data::kMCTrack, "mcreco" );
	scanner.ScanData( *pMCRecoTrack, (larlite::event_base*)ev_mctrack );

	auto pMCRecoShower = event.getValidHandle< std::vector<sim::MCShower> >( mcreco_track_tag );
	larlite::event_mcshower* ev_mcshower =
	  (larlite::event_mcshower*)llout.get_data( larlite::data::kMCShower, "mcreco" );
	scanner.ScanData( *pMCRecoShower, (larlite::event_base*)ev_mcshower );

	auto pSimCh = event.getValidHandle< std::vector<sim::SimChannel> >( simch_tag );
	larlite::event_simch* ev_simch =
	  (larlite::event_simch*)llout.get_data( larlite::data::kSimChannel, "simdrift" );
	scanner.ScanData( *pSimCh, (larlite::event_base*)ev_simch );
      }

      // opflash
      if ( HAS_OPFLASH ) {
// #ifdef ICARUS
//       auto popflashE = event.getValidHandle< std::vector<recob::OpFlash> >( opflashE_tag );
//       auto popflashW = event.getValidHandle< std::vector<recob::OpFlash> >( opflashW_tag );
//       larlite::event_opflash* ev_opflashE
// 	= (larlite::event_opflash*)llout.get_data( larlite::data::kOpFlash, "opflashCryoE" );
//       scanner.ScanData( *popflashE, ev_opflashE );
//       larlite::event_opflash* ev_opflashW
// 	= (larlite::event_opflash*)llout.get_data( larlite::data::kOpFlash, "opflashCryoW" );
//       scanner.ScanData( *popflashW, ev_opflashW );
// #endif
      }

      // CRT
      if ( HAS_CRT ) {
// #ifdef ICARUS
// 	auto pcrthit = event.getValidHandle< std::vector<sbn::crt::CRTHit> >( crthit_tag );
// 	larlite::event_crthit* ev_crthit
// 	  = (larlite::event_crthit*)llout.get_data( larlite::data::kCRTHit, "crthit" );
// 	scanner.ScanData( *pcrthit, ev_crthit );

// 	auto pcrttrack = event.getValidHandle< std::vector<sbn::crt::CRTTrack> >( crttrack_tag );
// 	larlite::event_crttrack* ev_crttrack
// 	  = (larlite::event_crttrack*)llout.get_data( larlite::data::kCRTTrack, "crttrack" );
// 	scanner.ScanData( *pcrttrack, ev_crttrack );
// #endif
      }
      
      llout.set_id( 0, 0, ientry );
      llout.next_event();
      // define a meta
      // we transfer data to image
      // we write output
      // if (true)
      // 	break;

      if ( ientry>=10 )
	break;
    }//end of event loop

  
  //geom->Print(std::cout);
  iolcv.finalize();
  llout.close();
    
  return 0;
} // main()

#endif // !__CLING__

