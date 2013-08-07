
#if !defined( __POINT_PROCESS_CORE_GAUSSIAN_POINT_PROCESS_UTILS_HPP__ )
#define __POINT_PROCESS_CORE_GAUSSIAN_POINT_PROCESS_UTILS_HPP__

#include <math-core/math_function.hpp>
#include <math-core/matrix.hpp>
#include <probability-core/distribution_utils.hpp>
#include <gsl/gsl_sf_erf.h>


namespace point_process_core {

  using namespace math_core;
  using namespace probability_core;


  // Description:
  // The likelihood of a negtaive observation (a region) given
  // a covariance ( so the mean is the function input / domain )
  class negative_observation_likelihood_for_mean_t
    : public math_function_t<nd_point_t,double>
  {
  public:
    double num_points_lambda;
    nd_aabox_t region;
    dense_matrix_t covariance;
    negative_observation_likelihood_for_mean_t
    ( const nd_aabox_t& region,
      const dense_matrix_t& covariance,
      const double& num_points_lambda )
      : region(region),
	covariance( covariance ),
	num_points_lambda( num_points_lambda )
    {}
    
    virtual
    double operator() ( const nd_point_t& mu ) const
    {
      // treat each dimension of mean independently
      double p = 1;
      for( size_t i = 0; i < mu.n; ++i ) {
	double x = mu.coordinate[i];
	double a = region.start.coordinate[i];
	double b = region.end.coordinate[i];
	double sig = to_eigen_mat(covariance)(i,i);
	double amass = gsl_sf_erfc( - ( a - x ) / ( sqrt(2.0) * sig ) );
	double bmass = gsl_sf_erfc( - ( b - x ) / ( sqrt(2.0) * sig ) );
	double diff = amass - bmass;
	double lik = exp( 0.5 * num_points_lambda * diff );
	p *= lik;
      }
      return p;
    }
  };


  // Descripiton:
  // The posterior distribution of a cluster  mean given both adata points
  // and negative observations (no longer conjugate hence we 
  // need the explicit posterior function )
  class gaussian_mixture_mean_posterior_t
    : public math_function_t<nd_point_t,double>
  {
  public:
    std::vector<nd_point_t> points;
    std::vector<nd_aabox_t> negative_observations;
    dense_matrix_t covariance;
    gaussian_distribution_t prior;
    poisson_distribution_t num_distribution;
    gaussian_distribution_t posterior_for_points_only;
    boost::shared_ptr<math_function_t<nd_point_t,double> > posterior;

    double scale;
    
    gaussian_mixture_mean_posterior_t
    ( const std::vector<nd_point_t>& points,
      const std::vector<nd_aabox_t>& negative_observations,
      const dense_matrix_t& cov,
      const poisson_distribution_t& num_distribution,
      const gaussian_distribution_t& prior )
      : points(points),
	negative_observations(negative_observations),
	covariance( cov ),
	num_distribution( num_distribution ),
	prior( prior )
    {
      calculate_posterior();
    }

    void calculate_points_only_posterior()
    {
      // the dimension of the poitns
      std::size_t dim = 0;
      if( points.empty() == false )
	dim = points[0].n;
      
      // invert the covariance to get a precision
      // (Not sure if this actually works for anything greater than 1D poitns!!)
      Eigen::MatrixXd prec_mat = to_eigen_mat( covariance ).inverse();
      
      // calculate the sum of the data points
      Eigen::VectorXd sum_vec( dim );
      for( size_t i = 0; i < points.size(); ++i ) {
	for( size_t k = 0; k < dim; ++k ) {
	  sum_vec(k) += points[i].coordinate[k];
	}
      }
      
      
      // Invert the prior's covariance to get a prior precision
      // and get the prior mean as Eign vector
      Eigen::MatrixXd prior_prec = to_eigen_mat( prior.covariance ).inverse();
      Eigen::VectorXd prior_mean = to_eigen_mat( prior.means );
      
      // Ok, now calculate the distribution over the new mixture mean
      Eigen::MatrixXd new_dist_prec = prec_mat * points.size() + prior_prec;
      Eigen::MatrixXd new_dist_cov = new_dist_prec.inverse();
      Eigen::VectorXd new_dist_mean 
	= new_dist_cov * ( prec_mat * sum_vec + prior_prec * prior_mean );
      
      // This is the posterior if you only include the observation points
      // and do NOT use the negative observations (so conjugate hence gaussian)
      gaussian_distribution_t new_dist;
      new_dist.dimension = dim;
      new_dist.means = to_vector( new_dist_mean ).component;
      new_dist.covariance = to_dense_mat( new_dist_cov );

      posterior_for_points_only = new_dist;
    }

    void calculate_posterior()
    {

      // cacluate the posterior using only the points and nopt
      // the negative regions
      calculate_points_only_posterior();

      // create the posterior math function 
      // (we will start with the point posterior and multiply by the likelihood
      // of the negaztive observatiosn )
      posterior = functions::gaussian_pdf( posterior_for_points_only );

      // Now, we need to multiply by the probability of *each* negative 
      // observation
      for( size_t i = 0; i < negative_observations.size(); ++i ) {
	boost::shared_ptr<math_function_t<nd_point_t,double> > neg_lik( new negative_observation_likelihood_for_mean_t( negative_observations[i], covariance, num_distribution.lambda ) );
	posterior = neg_lik * posterior;
      }

      // set the scale to the mean with only the points
      this->scale = pdf( point(posterior_for_points_only.means),
			 posterior_for_points_only );
    }

    double operator() ( const nd_point_t& mu ) const
    {
      return (*posterior)(mu);
    }
    
  };


}

#endif

