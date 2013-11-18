
#if !defined( __P2L_POINT_PROCESS_CORE_histogram_HPP__ )
#define __P2L_POINT_PROCESS_CORE_histogram_HPP__


#include "marked_grid.hpp"
#include <algorithm>

namespace point_process_core {


  //-------------------------------------------------------------------------
  

  // Description:
  // A histrogram stores a number of bins of a window of space and
  // a count of the number of elements seen within each bin.
  // This class builkds upon the marked grid to easily increment and 
  // decrement bin counts.
  template< class T >
  class histogram_t : public marked_grid_t<T>
  {
  public:
    
    // Descripotion:
    // Creates a new histogram with a set number of bins per dimension
    // over the given window
    histogram_t( const nd_aabox_t& window,
		 const size_t& bins_per_dim )
      : _bins_per_dim( bins_per_dim )
    {
      
      assert( !math_core::undefined(window) );
      
      // create the size (resolutions) for each of the dimansion of the
      // marked grid from the window sizes and num bins
      std::vector<double> resolutions;
      for( int i = 0; i < window.n; ++i ) {
	resolutions.push_back( (window.end.coordinate[i] - window.start.coordinate[i]) / (double)bins_per_dim );
      }
      
      // initialize the marked grid
      this->_init( window, window.start, resolutions );
    }
    
    
    // Description:
    // Increment the count for a given grid cell
    void increment_bin( const marked_grid_cell_t& cell, const T& inc = T(1) ) 
    {
      boost::optional<T> mark = this->operator()( cell );
      if( mark ) {
	this->set( cell, *mark + inc );
      } else {
	this->set( cell, inc );
      }
    }
    
    // Description:
    // increment count for the cell for the given point
    void increment_bin( const math_core::nd_point_t& p, const T& inc = T(1) )
    {
      this->increment_bin( this->cell( p ), inc );
    }

    // Description:
    // decrement the count for a given grid cell
    void decrement_bin( const marked_grid_cell_t& cell, const T& dec = T(1) ) 
    {
      boost::optional<T> mark = this->operator()( cell );
      if( mark ) {
	this->set( cell, *mark - dec );
      } else {
	this->set( cell, dec );
      }
    }

    // Description:
    // decrement count for the cell for the given point
    void decrement_bin( const math_core::nd_point_t& p, const T& dec = T(1) )
    {
      this->decrement_bin( this->cell( p ), dec );
    }

    // Description:
    // Retursn the total count in all bins in this histogram
    T total_count() const
    {
      T sum = T(0);
      for( auto cell : this->all_marked_cells() ) {
	sum += ( *this->operator()( cell ) );
      }
      return sum;
    }

    // Description:
    // Returns the number of bins per dimension
    size_t bins_per_dimension() const {
      return _bins_per_dim;
    }

  protected:

    // Description:
    // THe number of bins per dimension
    size_t _bins_per_dim;
    
  };
  
  
  //-------------------------------------------------------------------------


  // Description:
  // Creates a histogram from a set of data points
  template<typename T = double>
  histogram_t<T>
  create_histogram( const size_t& num_bins_per_dim,
		    const std::vector<math_core::nd_point_t>& samples )
  {
    histogram_t<T> hist( math_core::smallest_enclosing_box(samples),
			 num_bins_per_dim );
    for( auto s : samples ) {
      hist.increment_bin( s );
    }
    return hist;
  }

  //-------------------------------------------------------------------------
  

  // Description:
  // Normalize the histogram to contain probability masses rather than counts
  template<typename T, typename T_Result = double>
  histogram_t<T_Result>
  normalize_histogram( const histogram_t<T>& hist )
  {
    histogram_t<T_Result> norm( hist.window(), hist.bins_per_dimension() );

    // compute the total sum of counts
    T_Result sum = T_Result(0.0);
    for( auto cell : hist.all_marked_cells() ) {
      sum += ( *hist(cell) );
    }
    
    // create normalized histogram
    for( auto cell : hist.all_marked_cells() ) {
      norm.set( cell, T_Result( *hist(cell) / sum ) );
    }

    return norm;
  }


  //-------------------------------------------------------------------------
  
  
  // Description:
  // Computes the KL Divergence (Kullback-Leibler Divergence) between two
  // histograms when treated like distributions.
  // We use "absolute discounting" to nesure that the support of both
  // histograms is the same.  The epsilon term is added whenever one
  // histogram has a zero element. Then, the number of epsilons added is
  // removed from the "normalized" histogram values to get hte final
  // normalized histograms that are used for the kl-divergence
  template< typename TP, typename TQ, typename T_Result = double >
  T_Result kl_divergenge( const histogram_t<TP>& p_hist,
			  const histogram_t<TQ>& q_hist,
			  const T_Result& epsilon = T_Result(0.00001) )
  {
    // first, compute "pre" normalized histograms without doing any
    // smoothing of zero values
    histogram_t<T_Result> p_prenorm = normalize_histogram( p_hist );
    histogram_t<T_Result> q_prenorm = normalize_histogram( q_hist );

    // compute the set of elements with non-zero probability
    // from both histograms
    std::vector<marked_grid_cell_t> nonzero_cells = p_prenorm.all_marked_cells();
    std::vector<marked_grid_cell_t> q_cells = q_prenorm.all_marked_cells();
    nonzero_cells.insert( nonzero_cells.end(), q_cells.begin(), q_cells.end() );
    std::sort( nonzero_cells.begin(), nonzero_cells.end() );
    auto iter = std::unique( nonzero_cells.begin(),
			     nonzero_cells.end() );
    nonzero_cells.resize( std::distance( nonzero_cells.begin(), iter ) );

    // now smooth the "pre" normlaized histograms with the given epsilon
    // values
    histogram_t<T_Result> p = p_prenorm;
    histogram_t<T_Result> q = q_prenorm;
    size_t p_eps_count = 0;
    size_t q_eps_count = 0;
    std::vector<marked_grid_cell_t> p_nonzero_presmooth = p.all_marked_cells();
    std::vector<marked_grid_cell_t> q_nonzero_presmooth = q.all_marked_cells();
    
    // to smooth, first add epsilons for zero values, keep track of how many
    for( auto cell : nonzero_cells ) {
      boost::optional<T_Result> p_prob = p( cell );
      boost::optional<T_Result> q_prob = q( cell );
      
      // if we find a zero value, set it to epsilon and increment
      // the counts of epsilons used for p or q
      if( !p_prob ) {
	p_eps_count++;
	p.set( cell, epsilon );
      }
      if( !q_prob ) {
	q_eps_count++;
	q.set( cell, epsilon );
      }

      // one of the histograms should have non-zero, assert this
      assert( !( !p_prob && !q_prob ) ); 
    }

    // finish smoothing by subtracting the number of epsilons 
    // added from the other values in the normalized histogram
    // so that it is still normalized after the additional epsilons
    T_Result p_subnorm = ( p_eps_count * epsilon ) / p_nonzero_presmooth.size();
    T_Result q_subnorm = ( q_eps_count * epsilon ) / q_nonzero_presmooth.size();
    for( auto cell : p_nonzero_presmooth ) {
      p.decrement_bin( cell, p_subnorm );
    }
    for( auto cell : q_nonzero_presmooth ) {
      q.decrement_bin( cell, q_subnorm );
    }
    
    // now compute kl divergence
    T_Result kl = T_Result(0.0);
    for( auto cell : nonzero_cells ) {
      boost::optional<T_Result> p_prob = p( cell );
      boost::optional<T_Result> q_prob = q( cell );
      
      // we should have BOTH be non-zero because of smoothing
      assert( p_prob && q_prob );

      // compute ln( p/q ) * p
      kl += (*p_prob) * log( (*p_prob) / (*q_prob) );
    }

    // returns the kl divergence
    return kl;
  }
  

  //-------------------------------------------------------------------------
  

}


#endif

