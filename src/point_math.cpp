
#include "point_math.hpp"
#include <math-core/geom.hpp>


namespace point_process_core {


  //=======================================================================

  nd_point_t mean( const std::vector<nd_point_t>& points ) 
  {
    nd_point_t m;
    m.n = points[0].n;
    m.coordinate = std::vector<double>( m.n );
    
    // just take average along each coordinate
    for( size_t i = 0; i < points.size(); ++i ) {
      for( size_t dim = 0; (long)dim < m.n; ++dim ) {
	m.coordinate[dim] += points[i].coordinate[dim];
      }
    }
    for( size_t dim = 0; (long)dim < m.n; ++dim ) {
      m.coordinate[dim] /= points.size();
    }
    
    return m;
  }

  //=======================================================================

  double variance( const std::vector<nd_point_t>& points )
  {
    nd_point_t mean_p = mean( points );
    double sum = 0.0;
    for( size_t i = 0; i < points.size(); ++i ) {
      double diff = magnitude( points[i] - mean_p );
      sum += (diff * diff);
    }
    return sum / points.size();
  }

  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================


}
