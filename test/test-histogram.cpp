
#define BOOST_TEST_MODULE histogram
#include <boost/test/included/unit_test.hpp>

#include <point-process-core/histogram.hpp>
#include <probability-core/distributions.hpp>
#include <math-core/matrix.hpp>
#include <math-core/geom.hpp>
#include <math-core/io.hpp>
#include <iostream>
#include <sstream>

using namespace math_core;
using namespace probability_core;
using namespace point_process_core;


BOOST_AUTO_TEST_SUITE( test_suite_histogram )


// a simple fixture which generates some samples from a unit gaussian
struct fixture_unit_gaussian_samples
{
  fixture_unit_gaussian_samples() 
  {
    unit_gaussian.dimension = 1;
    unit_gaussian.means.push_back( 0.0 );
    unit_gaussian.covariance = diagonal_matrix( point( 1.0 ) );
    for( int i = 0; i < 1000; ++i ) {
      nd_point_t s = sample_from( unit_gaussian );
      if( i < 10 ) {
	samples_10.push_back( s );
      }
      if( i < 100 ) {
	samples_100.push_back( s );
      }	
      samples_1000.push_back( s );
    }
  }
  ~fixture_unit_gaussian_samples()
  {
  }
  gaussian_distribution_t unit_gaussian;
  std::vector<nd_point_t> samples_10;
  std::vector<nd_point_t> samples_100;
  std::vector<nd_point_t> samples_1000;
};


BOOST_FIXTURE_TEST_CASE( histogram_creation, fixture_unit_gaussian_samples )
{

  histogram_t<size_t> hist_10 = create_histogram<size_t>( 100, samples_10 );
  
  // make sure we have exatly ten counts
  BOOST_CHECK_EQUAL( hist_10.total_count(), size_t(10) );
  
  // amke sure every sample is inside the window
  for( auto s : samples_10 ) {
    BOOST_CHECK( is_inside( s, hist_10.window() ) );
  }

}


BOOST_FIXTURE_TEST_CASE( histogram_kl, fixture_unit_gaussian_samples )
{
  histogram_t<size_t> hist_10 = create_histogram<size_t>( 100, samples_10 );
  histogram_t<size_t> hist_100 = create_histogram<size_t>( 100, samples_100 );
  histogram_t<size_t> hist_1000 = create_histogram<size_t>( 100, samples_1000 );
  
  double kl_10_100 = kl_divergenge( hist_10, hist_100 );
  double kl_10_1000 = kl_divergenge( hist_10, hist_1000 );
  double kl_100_1000 = kl_divergenge( hist_100, hist_1000 );
  
  BOOST_CHECK_GT( kl_10_1000, kl_100_1000 );
}


BOOST_AUTO_TEST_SUITE_END()
