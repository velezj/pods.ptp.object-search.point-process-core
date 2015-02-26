
#if !defined( __POINT_PROCESS_CORE_ENTROPY_HPP__ )
#define __POINT_PROCESS_CORE_ENTROPY_HPP__

#include <math-core/types.hpp>
#include "point_process.hpp"

namespace point_process_core {


  // Descrioption:
  // A point process sampler
  // Given a state, returns a point set
  typedef std::vector<math_core::nd_point_t> (*point_process_sampler_t) ( void* state );


  // Description:
  // Parameters for extimating the entropy
  struct entropy_estimator_parameters_t
  {
    size_t num_samples;
    size_t num_samples_to_skip;
    double histogram_grid_cell_size;
    
    entropy_estimator_parameters_t()
      : num_samples( 100 ),
	num_samples_to_skip( 0 ),
	histogram_grid_cell_size( 1.0 )
    {}
  };
  
  // Description:
  // Estimate the entropy of a given point process from samples
  double estimate_entropy_from_samples
  ( const entropy_estimator_parameters_t& params,
    const math_core::nd_aabox_t& window,
    point_process_sampler_t sampler,
    void* state );

  // Description:
  // Estimate the entropy of a point process using samples from it
  double estimate_entropy_from_samples
  ( const entropy_estimator_parameters_t& params,
    boost::shared_ptr<mcmc_point_process_t>& process );
    

}


#endif
