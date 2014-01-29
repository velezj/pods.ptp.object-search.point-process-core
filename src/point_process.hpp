
#if !defined( __POINT_PROCESS_CORE_POINT_PROCESS_HPP__ )
#define __POINT_PROCESS_CORE_POINT_PROCESS_HPP__

#include <lcmtypes/p2l_math_core.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include "histogram.hpp"

namespace point_process_core {

  
  // Description:
  // A base class for all point processes
  class mcmc_point_process_t 
  {
  public:

    // Description:
    // Clones a point process
    // (this is usually so that we can compute a future state)
    virtual
    boost::shared_ptr<mcmc_point_process_t>
    clone() const = 0;

    // Description:
    // Turns on mcmc tracing with the given directory as the
    // trace root.
    virtual
    void trace_mcmc( const std::string& trace_dir ) = 0;
    virtual
    void trace_mcmc_off( ) = 0;

    // Description:
    // Retrusn the window for this point process
    virtual
    math_core::nd_aabox_t window() const = 0;

    // Description:
    // Returns the obsdervations for this process
    virtual 
    std::vector<math_core::nd_point_t>
    observations() const = 0;

    // Description:
    // Returns a point set sample from this process
    virtual
    std::vector<math_core::nd_point_t>
    sample() const = 0;

    // Description:
    // Runs a single step of MCMC sampling
    virtual
    void single_mcmc_step() = 0;
    
    // Description:
    // Returns a sample and steps the process by one
    virtual
    std::vector<math_core::nd_point_t>
    sample_and_step()
    {
      std::vector<math_core::nd_point_t> s = this->sample();
      this->single_mcmc_step();
      return s;
    }
    
    // Description:
    // Update with a new set of observations
    virtual
    void add_observations( const std::vector<math_core::nd_point_t>& obs ) = 0;

    // Description:
    // Add negative observation region
    virtual
    void add_negative_observation( const math_core::nd_aabox_t& region ) = 0;

    // Descripton:
    // Runs mcmc a number of iterations
    virtual
    void mcmc( const std::size_t& iterations, bool tick = false )
    {
      for( std::size_t i = 0; i < iterations; ++i ) {
	this->single_mcmc_step();
	if( tick ) {
	  std::cout << "." << i << "/" << iterations << "/";
	  std::cout.flush();
	}   
      }
    }


    // Description:
    // Prints a "shallow" one line trace (no newlines) for this
    // model parameters and state.
    virtual
    void print_shallow_trace( std::ostream& out ) const = 0;

    // Descripiton:
    // Plots the state of this point process.
    // The plot will have the given title.
    // The plot id is returned.
    // This uses the plot-server API for plotting, hence the
    // resulting ID allows fetching of the plot data itself.
    virtual
    std::string
    plot( const std::string& title ) const 
    { std::cout << "!!! DEFAULT NON-PLOT plot() called!" << std::endl;
      return ""; }


    // Description:
    // Computes an estimate for the intensity funciton for this point process
    // wihtin the given window, gridded with the given number of grids.
    // Optionally, you can specify the number of mcmc steps between samples
    // as well as the number of samples to use to compute the expected number
    // of points per grid region.
    virtual
    histogram_t<double>
    intensity_estimate( const math_core::nd_aabox_t& window,
			const size_t bins_per_dimension,
			const size_t num_samples_for_estimate = 1000,
			const size_t num_mcmc_iterations_between_samples = 1 )
    {
      histogram_t<double> hist( window, bins_per_dimension );
      for(size_t sample_i = 0; sample_i < num_samples_for_estimate; ++sample_i){
	std::vector<math_core::nd_point_t> sample = this->sample();
	for( math_core::nd_point_t p : sample ) {
	  hist.increment_bin( p );
	}
	this->mcmc( num_mcmc_iterations_between_samples, false );
      }
      // normalize the counts by the numbr of samples to get
      // average intensity
      for( auto cell : hist.all_cells() ) {
	if( hist(cell) ) {
	  hist.set( cell, *hist(cell) / (double)num_samples_for_estimate );
	}
      }
      return hist;
    }

  protected:


  };

}

#endif

