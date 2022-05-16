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
  std::vector<int> plane_indices;
  std::vector<int> data;
  //std::vector< TVector3 > pos_v;
  int dim1 = 0;
  int dim2 = 0;
  intersectiondata.Branch( "cryostat", &cryostat, "cryostat/I");
  intersectiondata.Branch( "tpc",      &tpc,      "tpc/I");
  intersectiondata.Branch( "plane_indices",   &plane_indices );
  intersectiondata.Branch( "dim1", &dim1 );
  intersectiondata.Branch( "dim2", &dim2 );
  intersectiondata.Branch( "data", &data );
  intersectiondata.Branch( "matrix_v", &matrix_v );
  //intersectiondata.Branch( "pos_v", &pos_v );

  for (int icryo=0; icryo<(int)geom->Ncryostats(); icryo++) {
    for (int itpc=0; itpc<(int)geom->NTPCs(icryo); itpc++) {

      int nplanes = geom->Nplanes(itpc,icryo);
      cryostat = icryo;
      tpc = itpc;
      
      // for each tpc, we need to know the intersections between planes
      for (int ip1=0; ip1<nplanes; ip1++) {
	int nwires_p1 = geom->Nwires( ip1, itpc, icryo );
	dim1 = nwires_p1;
	for (int ip2=ip1+1; ip2<nplanes; ip2++) {
	  int nwires_p2 = geom->Nwires( ip2, itpc, icryo );
	  dim2 = nwires_p2;
	  std::cout << "Make intersections for cryo=" << icryo << " tpc=" << itpc << " planes (" << ip1 << "," << ip2 << ")" << std::endl;

	  int otherplane = 0;
	  if ( ip1==0 && ip2==1 )
	    otherplane = 2;
	  if ( ip1==1 && ip2==2 )
	    otherplane = 0;
	  if ( ip1==0 && ip2==2 )
	    otherplane = 1;
	  
	  TMatrixD intersect_matrix( nwires_p1, nwires_p2 );
	  data.resize( dim1*dim2, 0 );
	  //pos_v.clear();
	  //pos_v.resize( dim1*dim2 );

	  for (int iwire1=0; iwire1<nwires_p1; iwire1<iwire1++) {
	    int ch1 = geom->PlaneWireToChannel( iwire1, ip1, itpc, icryo );
	    for (int iwire2=0; iwire2<nwires_p2; iwire2<iwire2++) {	    
	      int ch2 = geom->PlaneWireToChannel( iwire2, ip2, itpc, icryo );
	      TVector3 intersection(0,0,0);
	      bool crosses = geom->ChannelsIntersect( ch1, ch2, intersection );

	      if ( crosses ) {
		int nwires_p3 = geom->Nwires(otherplane,itpc,icryo);

		double wireco = geom->WireCoordinate( intersection, otherplane, itpc, icryo );
		int otherwire = wireco;
		if ( otherwire>=0 && otherwire<nwires_p3 ) {
		  data[ iwire1*nwires_p2 + iwire2 ] = otherwire+1;
		  intersect_matrix[iwire1][iwire2] = otherwire+1;
		  //pos_v[ iwire1*nwires_p2 + iwire2 ] = intersection;
		}
		else {
		  data[ iwire1*nwires_p2 + iwire2 ] = 0;
		  intersect_matrix[iwire1][iwire2] = 0;
		}
	      }
	      else {
		data[ iwire1*nwires_p2 + iwire2 ] = 0;
		intersect_matrix[iwire1][iwire2] = 0;
	      }
	    }
	  }
	  //matrix_v.emplace_back( std::move(intersect_matrix) );
	  //matrix = intersect_matrix;
	  plane_indices = std::vector<int>( { ip1, ip2, otherplane } );
	  matrix_v.clear();
	  matrix_v.push_back( intersect_matrix );
	  intersectiondata.Fill();	  
	}
      }
      plane_indices.clear();
    }
  }

  fout.Write();
  
  return 0;
}
