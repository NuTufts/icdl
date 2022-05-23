#include <iostream>
#include <vector>
#include <map>
#include <ctime>
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/DetectorProperties.h"
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/larflow3dhit.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"

int main( int nargs, char** argv )
{

  std::cout << "Generate Wire-Overlap Tensors" << std::endl;

  auto geom = larlite::larutil::Geometry::GetME( larlite::geo::kICARUS );
  auto detprop = larutil::DetectorProperties::GetME();

  // LOAD THE INTERSECTION DATA
  TFile fmatrices("output_icarus_wireoverlap_matrices.root");
  TTree* intersectiondata = (TTree*)fmatrices.Get("intersectdata");
  int cryostat = 0;
  int tpc = 0;
  std::vector<int>* p_plane_indices = 0;    
  int dim1 = 0;
  int dim2 = 0;
  std::vector< int >* p_data = 0;
  intersectiondata->SetBranchAddress( "cryostat", &cryostat );
  intersectiondata->SetBranchAddress( "tpc",      &tpc );
  intersectiondata->SetBranchAddress( "dim1",     &dim1 );
  intersectiondata->SetBranchAddress( "dim2",     &dim2 );  
  intersectiondata->SetBranchAddress( "plane_indices",  &p_plane_indices );
  intersectiondata->SetBranchAddress( "data", &p_data );

  std::vector< TMatrixD > matrix_list_v;
  std::map< std::vector<int>, int > m_planeid_to_tree_entry;
  
  int nentries = intersectiondata->GetEntries();
  for (int i=0; i<nentries; i++) {
    intersectiondata->GetEntry(i);
    std::vector<int> index_v = { cryostat, tpc, p_plane_indices->at(0), p_plane_indices->at(1), p_plane_indices->at(2) };
    std::cout << "matrix[" << i << "] cryo=" << cryostat << " tpc=" << tpc << " p1=" << index_v[2] << " p2=" << index_v[3] << std::endl;
    TMatrixD mat( dim1, dim2 ); // (nrows, ncols)
    std::cout << "  dim1= " << dim1 << " dim2=" << dim2 << std::endl;
    for (int i=0; i<dim1; i++) {
      for (int j=0; j<dim2; j++) {
	mat[i][j] = (*p_data)[ i*dim2 + j ];
      }
    }
    matrix_list_v.push_back( mat );
    m_planeid_to_tree_entry[ index_v ] = (int)matrix_list_v.size()-1;
  }

  std::cout << "LOAD MATRICES" << std::endl;

  // LOAD THE IMAGE DATA
  larcv::IOManager iolcv( larcv::IOManager::kREAD );
  iolcv.add_in_file( "../out_larcv.root" );
  iolcv.initialize();

  larlite::storage_manager ioll( larlite::storage_manager::kWRITE );
  ioll.set_out_filename("out_spacepoints.root");
  ioll.open();

  std::vector< std::vector<int> > plane_combos_v;
  plane_combos_v.push_back( std::vector<int>( {0,1,2} ) );
  plane_combos_v.push_back( std::vector<int>( {0,2,1} ) );

  std::cout << "[START ENTRY LOOP]" << std::endl;
  std::cin.get();
  
  for (int ientry=0; ientry<iolcv.get_n_entries(); ientry++) {
    iolcv.read_entry(ientry);
    larcv::EventImage2D* ev_adc
      = (larcv::EventImage2D*)iolcv.get_data(larcv::kProductImage2D, "wire");
    auto const& adc_v = ev_adc->as_vector();

    larlite::event_larflow3dhit* ev_hit =
      (larlite::event_larflow3dhit*)ioll.get_data( larlite::data::kLArFlow3DHit, "spacepoints" );

    int start_plane = 6;
    int icryo    = start_plane/12; // 12 images per cryo
    int itpc     = (start_plane-12*icryo)/3; // 3 images per tpc
    std::cout << "start_plane[" << start_plane << "] icryo=" << icryo << " itpc=" << itpc << std::endl;

    std::vector< std::vector<const TMatrixD*> > intersect_vv;

    for (int ii=0; ii<(int)plane_combos_v.size(); ii++) {
      int ipl1 = plane_combos_v[ii][0];
      int ipl2 = plane_combos_v[ii][1];
      int ipl3 = plane_combos_v[ii][2];
      std::vector<int> plane_index1 = { icryo, itpc, ipl1, ipl2, ipl3 };
      auto it1 = m_planeid_to_tree_entry.find( plane_index1 );
      if ( it1==m_planeid_to_tree_entry.end() )
	throw std::runtime_error("bad plane index1");      
      const TMatrixD& intersect1 = matrix_list_v.at( it1->second );

      std::vector<const TMatrixD*> intersect_v;
      intersect_v.push_back( &intersect1 );
      intersect_vv.push_back( intersect_v );
    }
      
    std::cout << "ready to rock" << std::endl;

    auto const& start_img = adc_v.at(start_plane);
    auto const& meta = start_img.meta();
    std::cout << meta.dump() << std::endl;
    int nrows = meta.rows();
    int ncols = meta.cols();

    std::vector< std::vector<int> > triplet_v; // (col1,col2,col3,row)

    std::time_t start = std::clock();
    
    for (int irow=0; irow<nrows; irow++) {
      float tick = meta.pos_y(irow);      
      // loop through the plane colums and get columns above thresh
      //std::cout << "irow[" << irow << "] -----------------------" << std::endl;
      std::vector< std::vector<int> > plane_cols(3);
      for (int iplane=0; iplane<3; iplane++) {
     	int splaneid = start_plane+iplane;
     	auto const& img = adc_v.at(splaneid);
	//std::cout << "iplane[" << iplane << "] icryo=" << icryo << " itpc=" << itpc << " iplaneid=" << planeid << std::endl;
     	for (int icol=0; icol<(int)img.meta().cols(); icol++) {
     	  if ( img.pixel(irow,icol)>10.0 )
     	    plane_cols[iplane].push_back(icol);
     	}
	//std::cout << "  above thresh pixels: " << plane_cols[iplane].size() << std::endl;
      }
      
      for (int ii=0; ii<(int)plane_combos_v.size(); ii++) {
	int ipl1 = plane_combos_v[ii][0];
     	int ipl2 = plane_combos_v[ii][1];
	int ipl3 = plane_combos_v[ii][2];
     	const TMatrixD& intersect = *(intersect_vv[ii][0]);
     	//std::cout << ii << ": " << ipl1 << "," << ipl2 << std::endl;
	//std::cout << "matrix: " << intersect.GetNrows() << "x" << intersect.GetNcols() << " elemns=" << intersect.GetNoElements() << std::endl;
 	for (auto idx1 : plane_cols[ipl1] ) {
	  for (auto idx2: plane_cols[ipl2] ) {
	    //std::cout << "(" << idx1 << "," << idx2 << ")" << std::endl;
	    if ( intersect[idx1][idx2]> 0 ) {	    
	      // register doublet
	      //std::cout << " register doublet: " << idx1 << "," << idx2 << std::endl;

	      // other plane wire
	      int otherwire = intersect[idx1][idx2]-1;
	      if ( adc_v[start_plane+ipl3].pixel(irow,otherwire)>0 ) {
		std::vector<int> triplet(4,0);
		triplet[ipl1] = idx1;
		triplet[ipl2] = idx2;
		triplet[ipl3] = otherwire;
		triplet[3] = irow;
		triplet_v.push_back( triplet );

		larlite::larflow3dhit hit;
		hit.resize(3);
		TVector3 pos(0,0,0);
		int ch1 = geom->PlaneWireToChannel( idx1, ipl1, itpc, icryo );
		int ch2 = geom->PlaneWireToChannel( idx2, ipl2, itpc, icryo );
		bool crosses = geom->ChannelsIntersect( ch1, ch2, pos );
		for (int i=0; i<3; i++)
		  hit[i] = pos[i];
		hit[0] = detprop->ConvertTicksToX(tick,0,itpc,icryo);
		ev_hit->emplace_back( std::move(hit) );
	      }
	    }
	  }
	}
	//std::cout << "trips registered this loop: " << trips_registered << std::endl;
	//std::cout << "tot triplets: " << triplet_v.size() << std::endl;
      }//end of plane combo loop
	
      //std::cout << "[enter] to continue." << std::endl;
      //std::cin.get();
    }// end of row loop
    std::cout << "TPC tot triplets: " << triplet_v.size() << std::endl;
    std::time_t end = std::clock();
    std::cout << "time elapsed: " << (float)(end-start)/CLOCKS_PER_SEC << std::endl;

    ioll.set_id( 0, 0, ientry );
    ioll.next_event();
  }//end of entry loop

  // GENERATE THE SPACEPOINTS
  ioll.close();
  
  return 0;
}
