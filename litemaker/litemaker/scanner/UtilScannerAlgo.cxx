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

  void UtilScannerAlgo::SaveLArProperties(const fhicl::ParameterSet& config,
					  const ::geo::GeometryCore& geom )
  {

    std::cout << "UtilScannerAlgo::SaveLArProperties" << std::endl;

    TTree* _larp_tree = new TTree("LArProperties", "Saved LAr Properties by UtilScannerAAlgo");

    // LArProperties setup
    auto const _larp = testing::setupProvider<detinfo::LArPropertiesStandard>(config.get<fhicl::ParameterSet>("services.LArPropertiesService"));
  
    // DetectorClocks setup
    auto const detClock = testing::setupProvider<detinfo::DetectorClocksStandard>(config.get<fhicl::ParameterSet>("services.DetectorClocksService"));
    auto const clockData = detClock->DataForJob();
  
    // DetectorProperties setup
    auto const detProp = testing::setupProvider<detinfo::DetectorPropertiesStandard>(config.get<fhicl::ParameterSet>("services.DetectorPropertiesService"),
										     detinfo::DetectorPropertiesStandard::providers_type{&geom,
											 static_cast<detinfo::LArProperties const*>(_larp.get()),
											 // TODO type cast is required until issue #18001 is solved
											 static_cast<detinfo::DetectorClocks const*>(detClock.get())});
    auto const propData = detProp->DataFor( clockData );


    //--- Fill Variables ---//
    std::vector< Int_t >         fEfield_cryoid;
    std::vector< Int_t >         fEfield_tpcid;
    std::vector< Int_t >         fEfield_planeid;
    std::vector< Double_t >      fEfield_field;

    for (geo::CryostatGeo const& cryostat: geom.IterateCryostats()) {
      int cryoid = (int)cryostat.ID().Cryostat;
      
      const unsigned int nTPCs = cryostat.NTPC();
      for(unsigned int t = 0;  t < nTPCs; ++t) {
	int tpcid = t;
	const geo::TPCGeo& tpc = cryostat.TPC(t);
	const unsigned int nPlanes = tpc.Nplanes();
	for(unsigned int p = 0; p < nPlanes; ++p) {
	  int planeid = p;
	  
	  fEfield_cryoid.push_back(  cryoid );
	  fEfield_tpcid.push_back(   tpcid );
	  fEfield_planeid.push_back( planeid );
	  fEfield_field.push_back( propData.Efield(planeid) );
	}// plane loop
      }//tpc loop
    }//cryo loop
    
    Double_t                         fTemperature      = propData.Temperature();
    Double_t                         fElectronlifetime = propData.ElectronLifetime(); ///< microseconds
    Double_t                         fRadiationLength  = _larp->RadiationLength();  ///< g/cm^2

    Double_t                         fArgon39DecayRate = _larp->Argon39DecayRate(); ///<  decays per cm^3 per second

    // Following parameters are for use in Bethe-Bloch formula for dE/dx.
    auto larp_pset = config.get< fhicl::ParameterSet >("services.LArPropertiesService");
    auto detp_pset = config.get< fhicl::ParameterSet >("services.DetectorPropertiesService");
    Double_t fZ = larp_pset.get<double>("AtomicNumber");       ///< Ar atomic number
    Double_t fA = larp_pset.get<double>("AtomicMass");         ///< Ar atomic mass (g/mol)
    Double_t fI = larp_pset.get<double>("ExcitationEnergy");   ///< Ar mean excitation energy (eV)
    Double_t fSa= detp_pset.get<double>("SternheimerA");       ///< Sternheimer parameter a
    Double_t fSk= detp_pset.get<double>("SternheimerK");       ///< Sternheimer parameter k
    Double_t fSx0 = detp_pset.get<double>("SternheimerX0");    ///< Sternheimer parameter x0
    Double_t fSx1 = detp_pset.get<double>("SternheimerX1");    ///< Sternheimer parameter x1
    Double_t fScbar = detp_pset.get<double>("SternheimerCbar");///< Sternheimer parameter Cbar

    // Optical parameters for Dar 
    std::vector<Double_t> fFastScintEnergies = larp_pset.get< std::vector<double> >("FastScintEnergies");
    std::vector<Double_t> fFastScintSpectrum = larp_pset.get< std::vector<double> >("FastScintSpectrum");
    std::vector<Double_t> fSlowScintEnergies = larp_pset.get< std::vector<double> >("SlowScintEnergies");
    std::vector<Double_t> fSlowScintSpectrum = larp_pset.get< std::vector<double> >("SlowScintSpectrum");
    std::vector<Double_t> fAbsLengthEnergies = larp_pset.get< std::vector<double> >("AbsLengthEnergies");
    std::vector<Double_t> fAbsLengthSpectrum = larp_pset.get< std::vector<double> >("AbsLengthSpectrum");
    std::vector<Double_t> fRIndexEnergies    = larp_pset.get< std::vector<double> >("RIndexEnergies"   );
    std::vector<Double_t> fRIndexSpectrum    = larp_pset.get< std::vector<double> >("RIndexSpectrum"   );
    std::vector<Double_t> fRayleighEnergies  = larp_pset.get< std::vector<double> >("RayleighEnergies" );
    std::vector<Double_t> fRayleighSpectrum  = larp_pset.get< std::vector<double> >("RayleighSpectrum" );

    bool fScintByParticleType = larp_pset.get<bool>("ScintByParticleType");

    Double_t fProtonScintYield        = larp_pset.get<double>("ProtonScintYield"     );
    Double_t fProtonScintYieldRatio   = larp_pset.get<double>("ProtonScintYieldRatio");
    Double_t fMuonScintYield          = larp_pset.get<double>("MuonScintYield"       );
    Double_t fMuonScintYieldRatio     = larp_pset.get<double>("MuonScintYieldRatio"  );
    Double_t fPionScintYield          = larp_pset.get<double>("PionScintYield"       );
    Double_t fPionScintYieldRatio     = larp_pset.get<double>("PionScintYieldRatio"  );
    Double_t fKaonScintYield          = larp_pset.get<double>("KaonScintYield"       );
    Double_t fKaonScintYieldRatio     = larp_pset.get<double>("KaonScintYieldRatio"  );
    Double_t fElectronScintYield      = larp_pset.get<double>("ElectronScintYield"   );
    Double_t fElectronScintYieldRatio = larp_pset.get<double>("ElectronScintYieldRatio");
    Double_t fAlphaScintYield         = larp_pset.get<double>("AlphaScintYield"      );
    Double_t fAlphaScintYieldRatio    = larp_pset.get<double>("AlphaScintYieldRatio" );  

    Double_t fScintResolutionScale = larp_pset.get<double>("ScintResolutionScale");
    Double_t fScintFastTimeConst   = larp_pset.get<double>("ScintFastTimeConst"  );
    Double_t fScintSlowTimeConst   = larp_pset.get<double>("ScintSlowTimeConst"  );
    Double_t fScintBirksConstant   = larp_pset.get<double>("ScintBirksConstant"  );
    Double_t fScintYield           = larp_pset.get<double>("ScintYield"          );
    Double_t fScintYieldRatio      = larp_pset.get<double>("ScintYieldRatio"     );
  
    bool fEnableCerenkovLight = larp_pset.get<bool>("EnableCerenkovLight");

    std::vector<std::string> fReflectiveSurfaceNames = larp_pset.get<std::vector<std::string> >("ReflectiveSurfaceNames");
    std::vector<Double_t> fReflectiveSurfaceEnergies = larp_pset.get<std::vector<double> >("ReflectiveSurfaceEnergies");;
    std::vector<std::vector<Double_t> > fReflectiveSurfaceReflectances = larp_pset.get<std::vector<std::vector<double> > >("ReflectiveSurfaceReflectances");
    std::vector<std::vector<Double_t> > fReflectiveSurfaceDiffuseFractions = larp_pset.get<std::vector<std::vector<double> > >("ReflectiveSurfaceDiffuseFractions");

    //--- Set TTree Branches ---//
    _larp_tree->Branch("fEfield_cryoid", "std::vector<Int_t>",  &fEfield_cryoid);
    _larp_tree->Branch("fEfield_tpcid",  "std::vector<Int_t>",  &fEfield_tpcid);
    _larp_tree->Branch("fEfield_planeid","std::vector<Int_t>",  &fEfield_planeid);
    _larp_tree->Branch("fEfield_field","std::vector<Double_t>", &fEfield_field);
    _larp_tree->Branch("fTemperature",&fTemperature,"fTemperature/D");
    _larp_tree->Branch("fElectronlifetime",&fElectronlifetime,"fElectronlifetime/D");
    _larp_tree->Branch("fRadiationLength",&fRadiationLength,"fRadiationLength/D");

    _larp_tree->Branch("fArgon39DecayRate",&fArgon39DecayRate,"fArgon39DecayRate/D");

    _larp_tree->Branch("fZ",&fZ,"fZ/D");
    _larp_tree->Branch("fA",&fA,"fA/D");
    _larp_tree->Branch("fI",&fI,"fI/D");
    _larp_tree->Branch("fSa",&fSa,"fSa/D");
    _larp_tree->Branch("fSk",&fSk,"fSk/D");
    _larp_tree->Branch("fSx0",&fSx0,"fSx0/D");
    _larp_tree->Branch("fSx1",&fSx1,"fSx1/D");
    _larp_tree->Branch("fScbar",&fScbar,"fScbar/D");

    _larp_tree->Branch("fFastScintSpectrum","std::vector<Double_t>",&fFastScintSpectrum);
    _larp_tree->Branch("fFastScintEnergies","std::vector<Double_t>",&fFastScintEnergies);
    _larp_tree->Branch("fSlowScintSpectrum","std::vector<Double_t>",&fSlowScintSpectrum);
    _larp_tree->Branch("fSlowScintEnergies","std::vector<Double_t>",&fSlowScintEnergies);
    _larp_tree->Branch("fRIndexSpectrum","std::vector<Double_t>",&fRIndexSpectrum);
    _larp_tree->Branch("fRIndexEnergies","std::vector<Double_t>",&fRIndexEnergies);
    _larp_tree->Branch("fAbsLengthSpectrum","std::vector<Double_t>",&fAbsLengthSpectrum);
    _larp_tree->Branch("fAbsLengthEnergies","std::vector<Double_t>",&fAbsLengthEnergies);
    _larp_tree->Branch("fRayleighSpectrum","std::vector<Double_t>",&fRayleighSpectrum);
    _larp_tree->Branch("fRayleighEnergies","std::vector<Double_t>",&fRayleighEnergies);

    _larp_tree->Branch("fScintByParticleType",&fScintByParticleType,"fScintByParticleType/O");

    _larp_tree->Branch("fProtonScintYield",&fProtonScintYield,"fProtonScintYield/D");
    _larp_tree->Branch("fProtonScintYieldRatio",&fProtonScintYieldRatio,"fProtonScintYieldRatio/D");
    _larp_tree->Branch("fMuonScintYield",&fMuonScintYield,"fMuonScintYield/D");
    _larp_tree->Branch("fMuonScintYieldRatio",&fMuonScintYieldRatio,"fMuonScintYieldRatio/D");
    _larp_tree->Branch("fPionScintYield",&fPionScintYield,"fPionScintYield/D");
    _larp_tree->Branch("fPionScintYieldRatio",&fPionScintYieldRatio,"fPionScintYieldRatio/D");
    _larp_tree->Branch("fKaonScintYield",&fKaonScintYield,"fKaonScintYield/D");
    _larp_tree->Branch("fKaonScintYieldRatio",&fKaonScintYieldRatio,"fKaonScintYieldRatio/D");
    _larp_tree->Branch("fElectronScintYield",&fElectronScintYield,"fElectronScintYield/D");
    _larp_tree->Branch("fElectronScintYieldRatio",&fElectronScintYieldRatio,"fElectronScintYieldRatio/D");
    _larp_tree->Branch("fAlphaScintYield",&fAlphaScintYield,"fAlphaScintYield/D");
    _larp_tree->Branch("fAlphaScintYieldRatio",&fAlphaScintYieldRatio,"fAlphaScintYieldRatio/D");

    _larp_tree->Branch("fScintYield",&fScintYield,"fScintYield/D");
    _larp_tree->Branch("fScintResolutionScale",&fScintResolutionScale,"fScintResolutionScale/D");
    _larp_tree->Branch("fScintFastTimeConst",&fScintFastTimeConst,"fScintFastTimeConst/D");
    _larp_tree->Branch("fScintSlowTimeConst",&fScintSlowTimeConst,"fScintSlowTimeConst/D");
    _larp_tree->Branch("fScintYieldRatio",&fScintYieldRatio,"fScintYieldRatio/D");  
    _larp_tree->Branch("fScintBirksConstant",&fScintBirksConstant,"fScintBirksConstant/D");

    _larp_tree->Branch("fEnableCerenkovLight",&fEnableCerenkovLight,"fEnableCerenkovLight/O");  

    _larp_tree->Branch("fReflectiveSurfaceNames","std::vector<std::string>",&fReflectiveSurfaceNames);
    _larp_tree->Branch("fReflectiveSurfaceEnergies","std::vector<Double_t>",&fReflectiveSurfaceEnergies);
    _larp_tree->Branch("fReflectiveSurfaceReflectances","std::vector<std::vector<Double_t> >",&fReflectiveSurfaceReflectances);
    _larp_tree->Branch("fReflectiveSurfaceDiffuseFractions","std::vector<std::vector<Double_t> >",&fReflectiveSurfaceDiffuseFractions);

    _larp_tree->Fill();

  }

  
  
}
}
