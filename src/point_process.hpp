
#if !defined( __POINT_PROCESS_CORE_POINT_PROCESS_HPP__ )
#define __POINT_PROCESS_CORE_POINT_PROCESS_HPP__

#include <lcmtypes/math_core.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>

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
	  std::cout << ".";
	  std::cout.flush();
	}   
      }
    }

  protected:


  };

}

#endif
