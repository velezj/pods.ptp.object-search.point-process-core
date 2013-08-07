
#include "entropy.hpp"
#include "marked_grid.hpp"
#include <math-core/geom.hpp>
#include <math-core/io.hpp>
#include <algorithm>
#include <iterator>
#include <iostream>


using namespace math_core;

namespace point_process_core {


  //=========================================================================

  double estimate_entropy_from_samples
  ( const entropy_estimator_parameters_t& params,
    const math_core::nd_aabox_t& window,
    point_process_sampler_t sampler,
    void* state )
  {
    
    // an array of the marked grids of point samples seen
    std::vector<marked_grid_t<size_t> > grids;
    std::vector<size_t> grid_counts;

    // sample a number of times from the sampler,
    // keep counts of the seen grids
    for( size_t i = 0; i < params.num_samples; ++i ) {
      
      // sample a point set
      std::vector<nd_point_t> sample = 
	sampler( state );

      // debug
      // std::cout << "  sampled process (#" << sample.size() << ")" << std::endl;
      // for( size_t j = 0; j < sample.size(); ++j ) {
      // 	std::cout << "    s[" << j << "]: " << sample[j] << std::endl;
      // }

      // mark the grid according to point set
      marked_grid_t<size_t> sample_grid( window, params.histogram_grid_cell_size );
      for( size_t i = 0; i < sample.size(); ++i ) {
	
	// debug
	//std::cout << "  about to mark " << i << std::endl;
	//std::cout << std::flush;
	
	boost::optional<size_t> mark = sample_grid( sample[i] );
	if( mark ) {
	  sample_grid.set( sample[i], *mark + 1 );
	} else {
	  sample_grid.set( sample[i], 1 );
	}
	
	// debug
	//std::cout << "    marking " << sample[i] << std::endl;
	
      }

      // debug
      //std::cout << "  marked grid ";

      // now see if we already have this grid
      std::vector<marked_grid_t<size_t> >::iterator fiter
	= std::find( grids.begin(), grids.end(), sample_grid );
      if( fiter != grids.end() ) {
	grid_counts[ std::distance( grids.begin(), fiter ) ] += 1;
	
	// debug
	//std::cout << "SAME!" << std::endl;
	
      } else {
	grids.push_back( sample_grid );
	grid_counts.push_back( 1 );
	
	// debug
	//std::cout << "new" << std::endl;

      }

      // skip some sampels if needed
      for( size_t skip_i = 0; skip_i < params.num_samples_to_skip; ++skip_i ) {
	sampler( state );
      }
    }

    // debug
    //std::cout << "  counts: ";

    // compute empirical entropy of the grids sample
    double entropy = 0;
    for( int i = 0; i < grid_counts.size(); ++i ) {
      double p = grid_counts[i] / (double)params.num_samples;
      entropy += p * log(p);
      
      // debug
      //std::cout << grid_counts[i] << " ";
      
    }
    
    // debug
    //std::cout << std::endl;

    return -entropy;
  }


  //=========================================================================

  double estimate_entropy_from_samples
  ( const entropy_estimator_parameters_t& params,
    boost::shared_ptr<mcmc_point_process_t>& process )
  {
   // an array of the marked grids of point samples seen
    std::vector<marked_grid_t<size_t> > grids;
    std::vector<size_t> grid_counts;

    // sample a number of times from the sampler,
    // keep counts of the seen grids
    for( size_t i = 0; i < params.num_samples; ++i ) {
      
      // sample a point set and step it
      std::vector<nd_point_t> sample = process->sample_and_step();

      // mark the grid according to point set
      marked_grid_t<size_t> sample_grid( process->window(), 
					 params.histogram_grid_cell_size );
      for( size_t i = 0; i < sample.size(); ++i ) {
	boost::optional<size_t> mark = sample_grid( sample[i] );
	if( mark ) {
	  sample_grid.set( sample[i], *mark + 1 );
	} else {
	  sample_grid.set( sample[i], 1 );
	}	
      }

      // now see if we already have this grid
      std::vector<marked_grid_t<size_t> >::iterator fiter
	= std::find( grids.begin(), grids.end(), sample_grid );
      if( fiter != grids.end() ) {
	grid_counts[ std::distance( grids.begin(), fiter ) ] += 1;	
      } else {
	grids.push_back( sample_grid );
	grid_counts.push_back( 1 );
      }

      // skip some mcmc steps if wanted
      for( size_t skip_i = 0; skip_i < params.num_samples_to_skip; ++skip_i ) {
	process->single_mcmc_step();
      }
    }

    // compute empirical entropy of the grids sample
    double entropy = 0;
    for( int i = 0; i < grid_counts.size(); ++i ) {
      double p = grid_counts[i] / (double)params.num_samples;
      entropy += p * log(p);
    }
    
    return -entropy;
  }

  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================



}
