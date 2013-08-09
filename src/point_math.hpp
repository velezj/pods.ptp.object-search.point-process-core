
#if !defined( __POINT_MATH_HPP__ )
#define __POINT_MATH_HPP__

#include <lcmtypes/p2l_math_core.hpp>
#include <lcmtypes/p2l_probability_core.hpp>


namespace point_process_core {

  using namespace math_core;
  using namespace probability_core;

  // Derwscription:
  // Take the mean of a set of points.
  // Uses eucledian distance
  nd_point_t mean( const std::vector<nd_point_t>& points );
  
  // Description:
  // Returns the variance of a set of points.
  // Uses eucledian distances
  double variance( const std::vector<nd_point_t>& points );


}


#endif

