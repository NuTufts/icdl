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
#include "lardataobj/RecoBase/OpFlash.h"

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

  bool HAS_MC  = false;
  bool HAS_CRT = false;

  // read FHiCL configuration from a configuration file:
  fhicl::ParameterSet config = lar::standalone::ParseConfiguration(configFile);
  
  // // set up message facility (always picked from "services.message")
  lar::standalone::SetupMessageFacility(config, "servicesdemo");

  //geometry setup
  auto geom = lar::standalone::SetupGeometry<icarus::ICARUSChannelMapAlg>(config.get<fhicl::ParameterSet>("services.Geometry"));
  std::cout << "larsoft geometry loaded" << std::endl;
  
  larutil::LArUtilConfig::SetDetector( larlite::geo::kICARUS );
  auto llgeom = larlite::larutil::Geometry::GetME();
  std::cout << "larlite geometry loaded" << std::endl;

  // litemaker scanner: converts larsoft data products into larlite ones
  larlite::ScannerAlgo scanner; 
  scanner.Register( "largeant",  larlite::data::kMCParticle );
  scanner.Register( "generator", larlite::data::kMCTruth );
  scanner.Register( "generator", larlite::data::kMCFlux );
  scanner.Register( "opflashCryoE",  larlite::data::kOpFlash );
  scanner.Register( "opflashCryoW",  larlite::data::kOpFlash );
  scanner.Register( "crthit",        larlite::data::kCRTHit );
  scanner.Register( "crttrack",      larlite::data::kCRTTrack );

  larlite::storage_manager llout( larlite::storage_manager::kWRITE );
  llout.set_out_filename( "out_larlite.root" );
  llout.open();

  // get parameters needed for image making
  unsigned int max_wires_per_plane   = geom->MaxWires();
  unsigned int max_ticks_per_channel = 4096;
  std::cout << "max wires per plane: " << max_wires_per_plane << std::endl;  

  std::vector<std::string> allInputFiles
    //= { "../../testdata/mcc2021b_bnbnu/prodcorsika_bnb_genie_protononly_overburden_icarus_gen_filter_g4_detsim_24027129_667_reco1.root"};
    //= {"/icarus/data/users/tmw/first_mc_files/prodcorsika_bnb_genie_protononly_overburden_icarus_gen_filter_g4_detsim_24027129_667_reco1.root"};
    = {"/icarus/data/users/tmw/first_mc_files/missing_hits_examples/BNB_run7033_11c_neutrinos_BNB_20211129T012418-stage0.root"};

  // OUTPUT SETUP
  larcv::IOManager iolcv( larcv::IOManager::kWRITE );
  iolcv.set_out_file( "out_larcv.root" );
  iolcv.initialize();

  int numEvents(0);

  // data tags we want to extract
  //std::string wireproducer = "daqTPCROI";
  std::string wireproducer = "roifinder";
  std::string stage_name = "stage0";
  //std::string stage_name = "MCstage0";
  std::vector< art::InputTag > wire_tag_v
    = {   art::InputTag(wireproducer,"PHYSCRATEDATATPCEW",stage_name),
	  art::InputTag(wireproducer,"PHYSCRATEDATATPCEE",stage_name),
	  art::InputTag(wireproducer,"PHYSCRATEDATATPCWE",stage_name),
	  art::InputTag(wireproducer,"PHYSCRATEDATATPCWW",stage_name) };

  art::InputTag const mcparticle_tag("largeant","","G4");
  art::InputTag const mctruth_tag("generator","","GenBNBbkgr");
  art::InputTag const mcflux_tag("generator","","GenBNBbkgr");
  art::InputTag const gtruth_tag("generator","","GenBNBbkgr");
  art::InputTag const crttrack_tag("crttrack","",stage_name);
  art::InputTag const crthit_tag("crthit","",stage_name);
  art::InputTag const opflashE_tag("opflashCryoE","",stage_name);
  art::InputTag const opflashW_tag("opflashCryoW","",stage_name);
  
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

      auto ev_wireout = (larcv::EventImage2D*)iolcv.get_data(larcv::kProductImage2D,"wire");

      std::vector< larcv::Image2D > adc_v;
      for (int iplane=0; iplane<24; iplane++) {
	// 2 cryo * 4 tpc * 3 planes = 24 planes total
	larcv::ImageMeta meta( max_wires_per_plane, max_ticks_per_channel, 
			       (int)max_ticks_per_channel, (int)max_wires_per_plane,
			       0.0, 0.0, iplane );
	larcv::Image2D img(meta);
	img.paint(0.0);
	std::cout << "blank made: " << meta.dump() << std::endl;
	adc_v.emplace_back( std::move(img) );
      }

      // we get the wire data
      for ( auto& wire_tag : wire_tag_v ) {
	auto pWireVec = event.getValidHandle< std::vector<recob::Wire> >(wire_tag);
	int nwires = pWireVec->size();
	std::cout << " num wires in event: " << nwires << std::endl;


	int nmaxsignal = 0;
	for ( auto const& wire : *pWireVec )  {
	  std::vector<float>  signal = wire.Signal();
	  raw::ChannelID_t   channel = wire.Channel();
	  std::vector<geo::WireID> wire_v = geom->ChannelToWire(channel);
	  if ( signal.size()>nmaxsignal )
	    nmaxsignal = signal.size();
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
	    std::cout << "CH " << channel << ": no wires=" << wire_v.size() << " num ticks=" << signal.size() << std::endl;
	    //}
	    img.copy( 0, wid, signal, std::min((int)signal.size(),(int)max_ticks_per_channel) );
	  }
	  //std::cout << "CH " << channel << ": no wires=" << wire_v.size() << std::endl;
	}
	//std::cout << "max signal size: " << nmaxsignal << std::endl;
      }//end of wire-tag loop
	
      for (auto& img : adc_v ) {

	// compress image first
	img.compress( max_ticks_per_channel/4, max_wires_per_plane, larcv::Image2D::kAverage );

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
	auto pMCParticleVec = event.getValidHandle< std::vector<simb::MCParticle> >(mcparticle_tag);
	larlite::event_mcpart* ev_mcpart = 
	  (larlite::event_mcpart*)llout.get_data( larlite::data::kMCParticle, "largeant" );
	scanner.ScanData( *pMCParticleVec, ev_mcpart );

	auto pMCTruthVec = event.getValidHandle< std::vector<simb::MCTruth> >(mctruth_tag);
	larlite::event_mctruth* ev_mctruth = 
	  (larlite::event_mctruth*)llout.get_data( larlite::data::kMCTruth, "generator" );
	scanner.ScanData( *pMCTruthVec, ev_mctruth );

	auto pMCFluxVec = event.getValidHandle< std::vector<simb::MCFlux> >(mctruth_tag);
	larlite::event_mcflux* ev_mcflux = 
	  (larlite::event_mcflux*)llout.get_data( larlite::data::kMCFlux, "generator" );
	scanner.ScanData( *pMCFluxVec, ev_mcflux );

	auto pGTruthVec = event.getValidHandle< std::vector<simb::GTruth> >(mctruth_tag);
	larlite::event_gtruth* ev_gtruth = 
	  (larlite::event_gtruth*)llout.get_data( larlite::data::kGTruth, "generator" );
	scanner.ScanData( *pGTruthVec, ev_gtruth );
      }

      // opflash
      auto popflashE = event.getValidHandle< std::vector<recob::OpFlash> >( opflashE_tag );
      auto popflashW = event.getValidHandle< std::vector<recob::OpFlash> >( opflashW_tag );
      larlite::event_opflash* ev_opflashE
	= (larlite::event_opflash*)llout.get_data( larlite::data::kOpFlash, "opflashCryoE" );
      scanner.ScanData( *popflashE, ev_opflashE );
      larlite::event_opflash* ev_opflashW
	= (larlite::event_opflash*)llout.get_data( larlite::data::kOpFlash, "opflashCryoW" );
      scanner.ScanData( *popflashW, ev_opflashW );

      // CRT
      if ( HAS_CRT ) {
	auto pcrthit = event.getValidHandle< std::vector<sbn::crt::CRTHit> >( crthit_tag );
	larlite::event_crthit* ev_crthit
	  = (larlite::event_crthit*)llout.get_data( larlite::data::kCRTHit, "crthit" );
	scanner.ScanData( *pcrthit, ev_crthit );

	auto pcrttrack = event.getValidHandle< std::vector<sbn::crt::CRTTrack> >( crttrack_tag );
	larlite::event_crttrack* ev_crttrack
	  = (larlite::event_crttrack*)llout.get_data( larlite::data::kCRTTrack, "crttrack" );
	scanner.ScanData( *pcrttrack, ev_crttrack );
      }


      llout.set_id( 0, 0, ientry );
      llout.next_event();
      // define a meta
      // we transfer data to image
      // we write output
      // if (true)
      // 	break;
    }

    
  //geom->Print(std::cout);
  iolcv.finalize();
  llout.close();
    
  return 0;
} // main()

#endif // !__CLING__

