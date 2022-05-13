#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "larlite/LArUtil/Geometry.h"

int main( int nargs, char** argv )
{

  std::cout << "Generate Wire-Overlap Tensors" << std::endl;

  auto geom = larlite::larutil::Geometry::GetME( larlite::geo::kICARUS );

  TFile fout("output_icarus_wireoverlap_matrices.root", "RECREATE");
  TTree intersectiondata( "intersectdata", "Wire Intersection Matrices for ICARUS" );
  std::vector< TMatrixD > matrix_v;
  int cryostat = 0;
  int tpc = 0;
  std::vector< std::vector<int> > plane_indices;;
  int plane2 = 0;
  intersectiondata.Branch( "cryostat", &cryostat, "cryostat/I");
  intersectiondata.Branch( "tpc",      &tpc,      "tpc/I");
  intersectiondata.Branch( "plane_indices",   &plane_indices );
  intersectiondata.Branch( "matrix_v", &matrix_v );

  for (int icryo=0; icryo<(int)geom->Ncryostats(); icryo++) {
    for (int itpc=0; itpc<(int)geom->NTPCs(icryo); itpc++) {

      int nplanes = geom->Nplanes(itpc,icryo);

      // for each tpc, we need to know the intersections between planes
      for (int ip1=0; ip1<nplanes; ip1++) {
	int nwires_p1 = geom->Nwires( ip1, itpc, icryo );
	for (int ip2=ip1+1; ip2<nplanes; ip2++) {
	  int nwires_p2 = geom->Nwires( ip2, itpc, icryo );

	  std::cout << "Make intersections for cryo=" << icryo << " tpc=" << itpc << " planes (" << ip1 << "," << ip2 << ")" << std::endl;
	  
	  std::vector<int> intersect_vec( nwires_p1*nwires_p2 );
	  TMatrixD intersect_matrix( nwires_p1, nwires_p2 );

	  for (int iwire1=0; iwire1<nwires_p1; iwire1<iwire1++) {
	    int ch1 = geom->PlaneWireToChannel( iwire1, ip1, itpc, icryo );
	    for (int iwire2=0; iwire2<nwires_p2; iwire2<iwire2++) {	    
	      int ch2 = geom->PlaneWireToChannel( iwire2, ip2, itpc, icryo );
	      TVector3 intersection;
	      bool crosses = geom->ChannelsIntersect( ch1, ch2, intersection );

	      intersect_vec[ iwire1*nwires_p2 + iwire2 ] = (crosses) ? 1 : 0;
	      intersect_matrix[iwire1][iwire2] = (crosses) ? 1 : 0;
	      
	    }
	  }
	  matrix_v.emplace_back( std::move(intersect_matrix) );
	  std::vector<int> pindices = { ip1, ip2 };
	  plane_indices.push_back( pindices );
	}
      }

      cryostat = icryo;
      tpc = itpc;
      
      intersectiondata.Fill();
      matrix_v.clear();
      plane_indices.clear();
    }
  }

  fout.Write();
  
  return 0;
}
