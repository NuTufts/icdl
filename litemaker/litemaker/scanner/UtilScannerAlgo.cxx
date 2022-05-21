#include "UtilScannerAlgo.h"

// ART
// #include "art/Framework/Principal/Event.h"
// #include "art/Framework/Principal/Handle.h"
// #include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
//#include "art/Framework/Services/Registry/ServiceMacros.h"

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesStandardTestHelpers.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
// - DetectorClocks
#include "lardataalg/DetectorInfo/DetectorClocksStandardTestHelpers.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
// - LArProperties
#include "lardataalg/DetectorInfo/LArPropertiesStandardTestHelpers.h"
#include "lardataalg/DetectorInfo/LArPropertiesStandard.h"

// ROOT
#include "TTree.h"

namespace litemaker {
namespace scanner {

  void UtilScannerAlgo::SaveDetectorProperties(const fhicl::ParameterSet& config,
					       const ::geo::GeometryCore& geom )
  {    

    std::cout << "UtilScannerAlgo::SaveDetectorProperties" << std::endl;

    // auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    // auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);

    // LArProperties setup
    auto const larp = testing::setupProvider<detinfo::LArPropertiesStandard>(config.get<fhicl::ParameterSet>("services.LArPropertiesService"));
  
    // DetectorClocks setup
    auto const detClock = testing::setupProvider<detinfo::DetectorClocksStandard>(config.get<fhicl::ParameterSet>("services.DetectorClocksService"));
    auto const clockData = detClock->DataForJob();
  
    // DetectorProperties setup
    auto const detProp = testing::setupProvider<detinfo::DetectorPropertiesStandard>(config.get<fhicl::ParameterSet>("services.DetectorPropertiesService"),
										     detinfo::DetectorPropertiesStandard::providers_type{&geom,
											 static_cast<detinfo::LArProperties const*>(larp.get()),
											 // TODO type cast is required until issue #18001 is solved
											 static_cast<detinfo::DetectorClocks const*>(detClock.get())});
    auto const propData = detProp->DataFor( clockData );
    
    
    TTree* _detp_tree = new TTree("DetectorProperties","");
    
    // //--- Fill Variables ---//
    Double_t fSamplingRate  = detinfo::sampling_rate(clockData);         ///< in ns
    Int_t    fTriggerOffset = detinfo::trigger_offset(clockData);       ///< in # of clock ticks
    Double_t fElectronsToADC = propData.ElectronsToADC();       ///< conversion factor for # of ionization electrons to 1 ADC count
    UInt_t   fNumberTimeSamples = propData.NumberTimeSamples(); ///< number of clock ticks per event
    UInt_t   fReadOutWindowSize = propData.ReadOutWindowSize(); ///< number of clock ticks per readout window
    Double_t fTimeOffsetU = propData.TimeOffsetU();             ///< time offsets to convert spacepoint
    Double_t fTimeOffsetV = propData.TimeOffsetV();             ///< coordinates to hit times on each
    Double_t fTimeOffsetZ = propData.TimeOffsetZ();             ///< view
    Double_t fXTicksCoefficient = propData.GetXTicksCoefficient(); ///< Parameters for x<-->ticks

    std::vector<Int_t>    fXTicksOffsets_cryoid;
    std::vector<Int_t>    fXTicksOffsets_tpcid;
    std::vector<Int_t>    fXTicksOffsets_planeid;
    std::vector<Double_t> fXTicksOffsets_offset;
    
    for (geo::CryostatGeo const& cryostat: geom.IterateCryostats()) {
      int cryoid = (int)cryostat.ID().Cryostat;
      
      const unsigned int nTPCs = cryostat.NTPC();
      for(unsigned int t = 0;  t < nTPCs; ++t) {
	int tpcid = t;
	const geo::TPCGeo& tpc = cryostat.TPC(t);
	const unsigned int nPlanes = tpc.Nplanes();
	for(unsigned int p = 0; p < nPlanes; ++p) {
	  int planeid = p;
	  
	  fXTicksOffsets_cryoid.push_back(  cryoid );
	  fXTicksOffsets_tpcid.push_back(   tpcid );
	  fXTicksOffsets_planeid.push_back( planeid );
	  fXTicksOffsets_offset.push_back( propData.GetXTicksOffset(planeid,tpcid,cryoid) );
	}// plane loop
      }//tpc loop
    }//cryo loop
    
    
    // OLD
    // for(unsigned int i=0; i<fXTicksOffsets.size(); ++i)
    //   fXTicksOffsets[i] = propData.GetXTicksOffset(i,0,0);
    
    //--- Set TTree Branches ---//
    _detp_tree->Branch("fSamplingRate",&fSamplingRate,"fSamplingRate/D");
    _detp_tree->Branch("fTriggerOffset",&fTriggerOffset,"fTriggerOffset/I");
    _detp_tree->Branch("fElectronsToADC",&fElectronsToADC,"fElectronsToADC/D");
    _detp_tree->Branch("fNumberTimeSamples",&fNumberTimeSamples,"fNumberTimeSamples/i");
    _detp_tree->Branch("fReadOutWindowSize",&fReadOutWindowSize,"fReadOutWindowSize/i");
    _detp_tree->Branch("fTimeOffsetU",&fTimeOffsetU,"fTimeOffsetU/D");
    _detp_tree->Branch("fTimeOffsetV",&fTimeOffsetV,"fTimeOffsetV/D");
    _detp_tree->Branch("fTimeOffsetZ",&fTimeOffsetZ,"fTimeOffsetZ/D");
    _detp_tree->Branch("fXTicksCoefficient",&fXTicksCoefficient,"fXTicksCoefficient/D");
    _detp_tree->Branch("fXTicksOffsets_cryoid", &fXTicksOffsets_cryoid );
    _detp_tree->Branch("fXTicksOffsets_tpcid", &fXTicksOffsets_tpcid );
    _detp_tree->Branch("fXTicksOffsets_planeid", &fXTicksOffsets_planeid );
    _detp_tree->Branch("fXTicksOffsets_offset", &fXTicksOffsets_offset );
    
    _detp_tree->Fill();
  
    return;
  }
  
  
}
}
